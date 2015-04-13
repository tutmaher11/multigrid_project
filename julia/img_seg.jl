using Color,FixedPointNumbers,Images,ImageView

#cd("Desktop/APPM_6640/mg_project_jl") # just make sure you run this in the same directory as your test images

function coarsenAMG(A,M,G,gamma,theta)

  # set Abar to be a matrix containing only strong connections of A
  Abar = A
  Arsum = [ sum(A[i,:]) - A[i,i] for i=1:M ]
  for i=1:M
    for j=1:M
      if (i == j || A[i,j] < theta*Arsum[i])
        Abar[i,j] = 0
      end
    end
  end

  # lambda[i] denotes the number of nonzero entries in column i of Abar
  lambda = zeros(M,1)
  for i=1:M
    lambda[i] = countnz(Abar[:,i])
  end

  # T keeps track of which set each node has been assigned to; T[i] = 1 means the i-th node is a C-pt, 2 means it's an F-pt,
  # means it's unassigned
  T = zeros(M,1)

  # set salient nodes as C-pts
  for i=1:M
    if G[i] < gamma
      T[i] = 1
      lambda[i] = 0
    end
  end

  while countnz(T) < M

    # this loop finds the smallest index satisfying T[j] = 0 and lambda[j] = maximum(lambda)
    j=1
    while T[j] != 0 || lambda[j] < maximum(lambda)
      j += 1
    end

    T[j] = 1
    lambda[j] = 0

    # K is the set of unassigned nodes that are strongly influenced by node j
    K = []
    for k=1:M
      if T[k] == 0 && Abar[k,j] > 0
        K = [K; k]
      end
    end

    for k in K
      T[k] = 2
      lambda[k] = 0

      #H is the set of nodes that strongly influence k
      H = []
      for h=1:M
        if T[h] == 0 && Abar[k,h] > 0
          H = [H; h]
        end
      end

      for h in H
        lambda[h] += 1
      end

    end

  end

  # C is the vector of  C-pts for current grid level
  C = []
  for i=1:M
    if T[i] == 1
      C = [C;i]
    end
  end

  return C

end

function imageVCycle(l,M,I,A,S,V,G,gamma,theta,atil,rho,beta,d1,sigma)

  println("size of graph: ", M)

  # do not allow salient segments on level if the mesh isn't coarse enough yet; 100000 is effectively infinity
  if l <= sigma
    G[:] = 100000
  end

  C = coarsenAMG(A,M,G,gamma,theta)
  Mn = length(C)
  ln = l + 1

  println("size of coarsened graph: ",Mn)

  # if no coarsening was obtained, we're done
  if Mn == M
    U = speye(M,M) # means every current node is a segment
    return U
  end

  # set the interpolation matrix
  P = zeros(M,Mn)
  counter = 1

  for i=1:M
    Anrsum = sum([A[i,C[k]] for k=1:Mn])
    if i in C
      P[i,counter] = 1
      counter += 1
    else
      P[i,:] = A[i,C[:]] / Anrsum
    end
  end
  P = sparse(P)

  # set the column-scaled interpolation matrix
  Ptil = zeros(M,Mn)
  Pcsum = [sum(P[:,j]) for j=1:Mn]
  for i=1:M
    for j=1:Mn
      Ptil[i,j] = P[i,j] / Pcsum[j]
    end
  end

  Ptil = sparse(Ptil)

  # set the coarse-level intensity vector
  In = Ptil'*I

  # do this to change the data to a regular data type; don't really understand why this is necessary,
  # but the following lines cause an error if we don't
  I = I*speye(M,M)

  # set the coarse level intensity variance measure
  Snc = (Ptil')*(I.^2) - In.^2
  Snf = Ptil'*S
  Sn = [Snf Snc]

  # define the coarse-level coupling matrix
  An = P'*A*P

  # rescale using coarse-level intensity
  for i=1:Mn
    for j=1:Mn
      An[i,j] = An[i,j]*exp(-atil*abs(In[i] - In[j]))
    end
  end

  # if ln >= rho, rescale using multilevel variance
  if ln >= rho
    for i=1:Mn
      for j=1:Mn
        An[i,j] = An[i,j]*exp(-beta*norm(Sn[i,:] - Sn[j,:]))
      end
    end
  end

  # Wn is the coarse-level weighted area matrix
  Wn = An

  #Vn is the coarse-level area matrix
  Vn = P'*V*P

  #Ln is the coarse-level weighted boundary length matrix
  Ln = -An
  for i=1:Mn
    Ln[i,i] = -sum([Ln[i,k] for k=1:Mn]) + Ln[i,i]
  end

  # Bn is the coarse-level boundary length matrix
  Bn = -Vn
  for i=1:Mn
    Bn[i,i] = -sum([Bn[i,k] for k=1:Mn]) + Bn[i,i]
  end

  #Gn is the coarse-level saliency vector
  Gn = zeros(Mn,1)
  for i=1:Mn
    if G[C[i]] == 0
      Gn[i] = 0
    else
      if (Ln[i,i] / Bn[i,i]) / (Wn[i,i] / Vn[i,i]) > gamma
        Gn[i] = (Ln[i,i] / Bn[i,i]) / (Wn[i,i] / Vn[i,i])
      else
        Gn[i] = 0
      end
    end
  end

  #recursively coarsen the graph
  Un = imageVCycle(ln,Mn,In,An,Sn,Vn,Gn,gamma,theta,atil,rho,beta,d1,sigma)

  #find the current segmentation matrix from the coarse-level segmentation matrix
  U = P*Un

  for i=1:M
    for j=1:size(U)[2]
      if U[i,j] < d1
        U[i,j] = 0
      elseif U[i,j] > 1 - d1
        U[i,j] = 1
      end
    end
  end

  return U

end

function img_seg_Driver(img,alpha,atil,beta,theta,gamma,d1,sigma,rho)

  #NOTE: this function currently only accepts square images, i.e. images with n x n, rather than n x m, pixels

  #first, convert the image to gray-scale
  img_gs = convert(Image{Gray},img)

  l = 1

  M,pass = size(img_gs)
  M = M^2

  # set 1st-level intensity vector
  I = []
  for j=1:sqrt(M)
    I = [I; img_gs[:,j]]
  end

  # A is the level 1 coupling matrix; only assign weights if nodes i and j are horizontal or vertical neighbors

  A = zeros(M,M)
  for i=2:M
    if mod(i-1,sqrt(M)) != 0
      A[i-1,i] = exp(-alpha*abs(I[i-1] - I[i]))
      A[i,i-1] = A[i-1,i]
    end
  end
  for i=sqrt(M)+1:M
      A[i-sqrt(M),i] = exp(-alpha*abs(I[i-sqrt(M)] - I[i]))
      A[i,i-sqrt(M)] = A[i-sqrt(M),i]
  end

  A = sparse(A)

  # S is the variance matrix
  S = zeros(M,1)

  # L is the weighted boundary length matrix (graph Laplacian)
  L = -A
  Arsum = [sum(A[i,:]) - A[i,i] for i=1:M]
  for i=1:M
    L[i,i] = Arsum[i]
  end

  # V is the area matrix
  V = zeros(M,M)
  for i=2:M
    if mod(i-1,sqrt(M)) != 0
      V[i-1,i] = 1
      V[i,i-1] = 1
    end
  end
  for i=sqrt(M)+1:M
      V[i-sqrt(M),i] = 1
      V[i,i-sqrt(M)] = 1
  end

  V = sparse(V)

  # B is the boundary length matrix
  B = -V
  Vrsum = [sum(V[i,:]) - V[i,i] for i=1:M]
  for i=1:M
        B[i,i] = Vrsum[i]
  end

  # G is the saliency vector
  G = zeros(M,1)
  for i=1:M
    G[i] = L[i,i] / B[i,i]
  end

  U = imageVCycle(l,M,I,A,S,V,G,gamma,theta,atil,rho,beta,d1,sigma)

  # assign pixels uniquely to segments
  for i=1:M
    j = 1
    while U[i,j] != maximum(U[i,:]) # this assigns the pixel to the first segment which has the highest U[i,j] value
      j += 1
    end
    if j == 1
      U[i,1] = 1
      U[i,2:end] = 0
    elseif j == size(U)[2]
      U[i,end] = 1
      U[i,1:end-1] = 0
    else
      U[i,1:j-1] = 0
      U[i,j] = 1
      U[i,j+1:end] = 0
    end
  end

  return U

end

function img_seg_proc(img,U)

  # convert the image into a grayscale version that retains the data structure of a color image
  img_gs = convert(Image{RGB},convert(Image{Gray},img))

  #separate the image into a three-dimensional array; one matrix each for R, G, B
  img_gs = separate(img_gs)

  # convert each column of U into a matrix representing a grayscale image of the segment
  pix_num,seg_num = size(U)
  pic_width = convert(Int64,sqrt(pix_num))

  segs = zeros(pic_width,pic_width,seg_num)

  for k=1:seg_num
    for j=1:pic_width
      segs[:,j,k] = U[pic_width*(j-1)+1:pic_width*(j),k]
    end
    segs[:,:,k] = segs[:,:,k]' # picture coordinates of original image are given in opposite order when you convert
                               # a grayscale image to a color image? Or something.
  end

  for k=1:seg_num
    # convert each segment matrix into a blue-tinted image
    segs_img = zeros(pic_width,pic_width,3)
    segs_img[:,:,3] = 0.5*segs[:,:,k]
    segs_img = convert(Image{RGB},segs_img)
    segs_img = separate(segs_img)

    # add original picture and segment picture together and scale; result is a picture with just the segment tinted blue
    new_img = segs_img + img_gs
    new_img = new_img / maximum(new_img)

    view(new_img)

  end

end

img = imread("peppers_50.jpg")
U = img_seg_Driver(img,10,10,10,0.1,0.1,0.15,5,1)
img_seg_proc(img,U)
