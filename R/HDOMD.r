## Construction of high-dimensional orthogonal maximin distance designs with flexible numbers of runs and factors
## Author: Xu He (hexu@amss.ac.cn) and Fasheng Sun


HDOM2 <- function(n,p) {  # Construction of two-level designs

	BestSep = 0 
	n1list = c(2,seq(4,n/2,4))
	n1list = n1list[ floor(n / n1list) * n1list == n ]
	n2list = n / n1list
	
	H = Hadamard_Matrix(order=n)
	for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
	BestD = H[,(n-p+1):n,drop=FALSE]/2/2 + 1/2 
	BestSep = min(dist(BestD,method="manhattan"))
	Bestn1 = 1
	Bestn2 = n
	Besttildep = n
	Bestbarp = p
	Bestq = p
	if(BestSep<2) if(floor(log2(n))==log2(n)) if(p>log2(n)) {
		barP = c(n - 2^(1:(log2(n)-1)) , n)
		if(floor(log2(n)/2)*2==log2(n))  barP = c(barP,n-1)
		if(floor(log2(n)/2)*2!=log2(n))  barP = c(barP,2)
		E = H[, barP, drop=FALSE]
		Ec = H[, -barP, drop=FALSE]
		if(dim(E)[2]<p) E = cbind(Ec[,(dim(Ec)[2]-p+dim(E)[2]+1):dim(Ec)[2],drop=FALSE],E)
		BestD = E/2/2 + 1/2 
		BestSep = min(dist(BestD,method="manhattan"))
		Bestbarp = log2(n)+1
		Bestq = log2(n)+1
	}
	
	for(indexn in 1:length(n1list)){
		n1 = n1list[indexn]
		n2 = n2list[indexn] 
		if(n2>2) if(floor(n2/4)*4!=n2) next

		A = Hadamard_Matrix(order=n1)
		if(n1==1) A = matrix(1,1,1)
		if(length(A)!=n1^2) next
		for(i in 1:n1) if(A[i,1]==-1) A[i,] = -A[i,]
		A = A[,-1,drop=FALSE]
		
		n2b = n2*2
		while(floor(n2b/16)*16==n2b)  n2b = n2b/2
		if(n2b==8) n2b = 2
		if(n2b==4) n2b = 2
		H = Hadamard_Matrix(order=n2b/2)
		if(n2b>2) if(length(H)==1)  next 
		if(length(H)>1) for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
		while(length(H)!=n2^2)  H = rbind( cbind(H,H), cbind(H,-H) ) 
		
		for(i in 1:9) {
			if(i>=1) if(i<=5) {
				if(n2-(i-1)<=0) next
				B = H[,-(1:(i-1)),drop=FALSE]
				if(i==1) B = H
				Bc = H[,(1:(i-1)),drop=FALSE]
				if(i==1) Bc = matrix(0,dim(H)[1],0)
			}
			if(i>=6) if(i<=8) {
				if(n2/2-(i-6)<=0) next
				B = H[,-(1:(n2/2+i-6)),drop=FALSE]
				Bc = H[,1:(n2/2+i-6),drop=FALSE]
			}
			if(i==9) {
				if(log2(n2)!=floor(log2(n2))) next
				if(log2(n2)<=3) next
				barP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
				if(floor(log2(n2)/2)*2==log2(n2))  barP = c(barP,n2-1)
				if(floor(log2(n2)/2)*2!=log2(n2))  barP = c(barP,2)
				B = H[, barP, drop=FALSE]
				Bc = H[, -barP, drop=FALSE]
			}
			if(dim(B)[2]<floor((p-n2-1)/(n1-1))) next
			if(dim(B)[2]>ceiling(p/(n1-1))) next
			
			barp = p-(n1-1)*dim(B)[2]
			if(barp<0) barp = 0
			if(barp>n2-1) barp = n2-1 
			for(j in 1:2) {
				if(j==1) if(barp==0) {
					C = matrix(0,n2,0)
					Cc = H[,-1,drop=FALSE]
				}
				if(j==1) if(barp>0) {
					C = H[,(n2-barp+1):n2,drop=FALSE]
					Cc = H[,c(-1,-((n2-barp+1):n2)),drop=FALSE]
				}
				if(j==2) if(barp>0) {
					if(log2(n2)!=floor(log2(n2))) next
					if(log2(n2)<=3) next
					barP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
					if(floor(log2(n2)/2)*2==log2(n2))  barP = c(barP,n2-1)
					if(floor(log2(n2)/2)*2!=log2(n2))  barP = c(barP,2)
					C = H[, barP, drop=FALSE]
					Cc = H[, c(-1,-barP), drop=FALSE]
				}

				E = cbind( kronecker(A,B), kronecker(matrix(1,n1,1),C) )
				D = E/2/2 + 1/2 

				if(dim(Bc)[2]>0) Ec = cbind( kronecker(A,Bc), kronecker(matrix(1,n1,1),Cc) )
				if(dim(Bc)[2]==0) Ec = kronecker(matrix(1,n1,1),Cc)
				Dc = Ec/2/2 + 1/2
				if(dim(D)[2]>p) D = D[,(dim(D)[2]-p+1):dim(D)[2]]
				if(dim(D)[2]<p) { 
					D = cbind( Dc[ ,(dim(Ec)[2]-p+dim(D)[2]+1):dim(Ec)[2],drop=FALSE], D)
				}
				Sep = min(dist(D,method="manhattan"))
				FinalSep = Sep 
				
				if(FinalSep>BestSep) {
					BestSep = FinalSep
					BestD = D
					Bestn1 = n1
					Bestn2 = n2
					Besttildep = dim(B)[2]
					Bestbarp = dim(C)[2]
					Bestq = (n1-1)*dim(B)[2]+dim(C)[2]
				}
			}
		}
	}

	if(BestSep==-1) return(NULL)
	return(BestD)
}


R2 = rbind( c(2,1), c(-1,2) )/sqrt(5)

HDOM4 <- function(n,p) {  # Construction of four-level designs

	BestSep = 0 
	n1list = c(2,seq(4,n/2,4))
	n1list = n1list[ floor(n / n1list) * n1list == n ]
	n2list = n / n1list
	
	H = Hadamard_Matrix(order=n)
	for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
	E = H[,(n-ceiling(p/2)*2+1):n,drop=FALSE]
	W = E
	for(k in 1:(dim(W)[2]/2)) W[,(2*k-1):(2*k)] = E[,(2*k-1):(2*k)] %*% R2
	BestD = W[,(dim(W)[2]-p+1):dim(W)[2]]*sqrt(5)/4/2 + 1/2
	BestSep = min(dist(BestD))
	Bestn1 = 1
	Bestn2 = n
	Besttildep = n
	Bestbarp = p
	Bestq = p
	if(BestSep<2) if(floor(log2(n))==log2(n)) if(p>log2(n)) {
		barP = c(n - 2^(1:(log2(n)-1)) , n)
		if(floor(log2(n)/2)*2==log2(n))  barP = c(n-3,barP,n-1)
		if(floor(log2(n)/2)*2!=log2(n))  barP = c(barP,2)
		E = H[, barP, drop=FALSE]
		Ec = H[, -barP, drop=FALSE]
		if(dim(E)[2]<ceiling(p/2)*2) E = cbind(Ec[,(dim(Ec)[2]-ceiling(p/2)*2+dim(E)[2]+1):dim(Ec)[2],drop=FALSE],E)
		W = E
		for(k in 1:(dim(W)[2]/2)) W[,(2*k-1):(2*k)] = E[,(2*k-1):(2*k)] %*% R2
		D = W[,(dim(W)[2]-p+1):dim(W)[2]]*sqrt(5)/4/2 + 1/2
		Sep = min(dist(D))
		if(Sep>BestSep) {
			BestD = D
			BestSep = Sep
			Bestbarp = log2(n)+1
			Bestq = log2(n)+1
		}
	}
	
	for(indexn in 1:length(n1list)){
		n1 = n1list[indexn]
		n2 = n2list[indexn] 
		if(n2>2) if(floor(n2/4)*4!=n2) next

		A = Hadamard_Matrix(order=n1)
		if(n1==1) A = matrix(1,1,1)
		if(length(A)!=n1^2) next
		for(i in 1:n1) if(A[i,1]==-1) A[i,] = -A[i,]
		A = A[,-1,drop=FALSE]
		
		n2b = n2*2
		while(floor(n2b/16)*16==n2b)  n2b = n2b/2
		if(n2b==8) n2b = 2
		if(n2b==4) n2b = 2
		H = Hadamard_Matrix(order=n2b/2)
		if(n2b>2) if(length(H)==1)  next 
		if(length(H)>1) for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
		while(length(H)!=n2^2)  H = rbind( cbind(H,H), cbind(H,-H) ) 
		
		for(i in 1:9) {
			if(i>=1) if(i<=5) {
				if(n2-(i-1)<=0) next
				B = H[,-(1:(i-1)),drop=FALSE]
				if(i==1) B = H
				Bc = H[,(1:(i-1)),drop=FALSE]
				if(i==1) Bc = matrix(0,dim(H)[1],0)
			}
			if(i>=6) if(i<=8) {
				if(n2/2-(i-6)<=0) next
				B = H[,-(1:(n2/2+i-6)),drop=FALSE]
				Bc = H[,1:(n2/2+i-6),drop=FALSE]
			}
			if(i==9) {
				if(log2(n2)!=floor(log2(n2))) next
				if(log2(n2)<=3) next
				barP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
				if(floor(log2(n2)/2)*2==log2(n2))  barP = c(barP,n2-1)
				if(floor(log2(n2)/2)*2!=log2(n2))  barP = c(barP,2)
				B = H[, barP, drop=FALSE]
				Bc = H[, -barP, drop=FALSE]
			}
			if(dim(B)[2]<floor((p-n2-1)/(n1-1))) next
			if(dim(B)[2]>ceiling(p/(n1-1))) next
			
			barp = p-(n1-1)*dim(B)[2]
			if(barp<0) barp = 0
			if(barp>n2-1) barp = n2-1 
			for(j in 1:2) {
				if(j==1) if(barp==0) {
					C = matrix(0,n2,0)
					Cc = H[,-1,drop=FALSE]
				}
				if(j==1) if(barp>0) {
					C = H[,(n2-barp+1):n2,drop=FALSE]
					Cc = H[,c(-1,-((n2-barp+1):n2)),drop=FALSE]
				}
				if(j==2) if(barp>0) {
					if(log2(n2)!=floor(log2(n2))) next
					if(log2(n2)<=3) next
					barP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
					if(floor(log2(n2)/2)*2==log2(n2))  barP = c(barP,n2-1)
					if(floor(log2(n2)/2)*2!=log2(n2))  barP = c(barP,2)
					C = H[, barP, drop=FALSE]
					Cc = H[, c(-1,-barP), drop=FALSE]
				}

				E = cbind( kronecker(A,B), kronecker(matrix(1,n1,1),C) )
				if(dim(Bc)[2]>0) Ec = cbind( kronecker(A,Bc), kronecker(matrix(1,n1,1),Cc) )
				if(dim(Bc)[2]==0) Ec = kronecker(matrix(1,n1,1),Cc)
				qu = max(c(ceiling(p/2)*2,ceiling(dim(E)[2]/2)*2))
				if(dim(E)[2] < qu)  E = cbind(Ec[,(dim(Ec)[2]-qu+dim(E)[2]+1):dim(Ec)[2],drop=FALSE],E)
				W = E
				for(k in 1:(dim(W)[2]/2)) W[,(2*k-1):(2*k)] = E[,(2*k-1):(2*k)] %*% R2
				D = W[,(dim(W)[2]-p+1):dim(W)[2]]*sqrt(5)/4/2 + 1/2
				FinalSep = min(dist(D))
				
				if(FinalSep>BestSep) {
					BestSep = FinalSep
					BestD = D
					Bestn1 = n1
					Bestn2 = n2
					Besttildep = dim(B)[2]
					Bestbarp = dim(C)[2]
					Bestq = (n1-1)*dim(B)[2]+dim(C)[2]
				}
			}
		}
	}

	if(BestSep==-1) return(NULL)
	return(BestD)
}



R3 = rbind( c(4,-2,-1,0), c(2,4,0,1), c(1,0,4,-2), c(0,-1,2,4) )/sqrt(21)

R4 = rbind( c(8,-4,-2,1), c(4,8,-1,-2), c(2,-1,8,-4), c(1,2,4,8) )/sqrt(85)

gpK4 <- function(J,K) {
	if(length(K)==0) return(NULL)
	if(length(K)>4)  return( rbind( gpK4(J,K[1:4]), gpK4(J,K[5:length(K)]) ) )
	if(length(intersect(0,J))==0)  return( cbind(J,J,J,J[c(length(J),1:(length(J)-1))],K[1],K[2],K[3],K[4]) ) 
	if(length(intersect(0,J))>0) { 
		if(length(J)>4)  return( rbind( gpK4(J[1:(length(J)-2)],K), gpK4(J[(length(J)-1):length(J)],K) ) )
		if(length(J)==3) return( rbind( c(J[1],J[1],J[2],J[3],K[1],K[2],K[1],K[1]), 
			c(J[1],J[1],J[2],J[3],K[3],K[4],K[3],K[4]), c(J[2],J[2],J[3],J[3],K[2],K[4],K[2],K[3]) ) )
		if(length(J)==4) return( cbind(J[1],J[2],J[3],J[4],K,K,K,K[c(length(K),1:(length(K)-1))]) )
	}
}

gpJ4 <- function(J,K)  { if(length(J)==0) return(NULL); TTT = gpK4(K-1,J+1); return( cbind(TTT[,5:8]-1,TTT[,1:4]+1) ); }

# gpK2 <- function(J,K) {
	# if(length(J)>6)  return( rbind( gpK2(J[1:(length(J)-4)],K), gpK2(J[(length(J)-3):length(J),K) ) )
	# if(length(J)==4) return( cbind(J[1],J[2],J[3],J[4],K,K,K,K[2:1]) )
	# if(length(J)==6) return( rbind( c(J[1],J[1],J[2],J[3],K[1],K[2],K[1],K[1]), 
		# c(J[4],J[4],J[5],J[6],K[1],K[2],K[1],K[2]), c(J[2],J[3],J[5],J[6],K[2],K[2],K[2],K[1]) ) )
# }


HDOM16 <- function(n,p) {  # Construction of sixteen-level designs
	if(ceiling(p/4)*4>=n) return(NULL)
	BestSep = -1 
	n1list = c(seq(4,n/4,4))
	n1list = n1list[ floor(n / n1list) * n1list == n ]
	n2list = n / n1list
		
	for(indexn in 1:length(n1list)){
		n1 = n1list[indexn]
		n2 = n2list[indexn] 
		if(n2>2) if(floor(n2/4)*4!=n2) next

		A = Hadamard_Matrix(order=n1)
		if(n1==1) A = matrix(1,1,1)
		if(length(A)!=n1^2) next
		for(i in 1:n1) if(A[i,1]==-1) A[i,] = -A[i,]
		
		n2b = n2*2
		while(floor(n2b/16)*16==n2b)  n2b = n2b/2
		if(n2b==8) n2b = 2
		if(n2b==4) n2b = 2
		H = Hadamard_Matrix(order=n2b/2)
		if(n2b>2) if(length(H)==1)  next 
		if(length(H)>1) for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
		while(length(H)!=n2^2)  H = rbind( cbind(H,H), cbind(H,-H) ) 
		Efull = kronecker(A,H)
		
		for(tildep in floor((p-n2+1)/(n1-1)) : ceiling(p/(n1-1)) )  for(j2 in 1:2)  {
			if(tildep < 2) next
			if(tildep > n2) next
			
			j = 2
			if(log2(n2)!=floor(log2(n2))) j = 1
			if(log2(n2)+1>tildep) j = 1
			if(tildep>=n2/4+2) j = 1
			if(j==1)  tildeP = (n2-tildep+1):n2
			if(j==2) {
				tildeP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
				if(floor(log2(n2)/2)*2==log2(n2))  tildeP = c(tildeP,n2-1)
				if(floor(log2(n2)/2)*2!=log2(n2))  tildeP = c(tildeP,2)
				if(length(tildeP)<tildep) tildeP = c(tildeP,setdiff(1:n2,tildeP)[(n2-tildep+1):(n2-length(tildeP))])
				tildeP = sort(tildeP) 
			}
			
			barp = p - (n1-1)*tildep
			if(tildep>=n2-2)  barp = ceiling(max(c(p,(n1-1)*tildep))/4)*4 -(n1-1)*tildep;
			if(barp>min(c(n2-1,(n1-3)*tildep+2)))  barp = min(c(n2-1,(n1-3)*tildep+2));
			if(barp<0)  barp = 0;
			addp = ceiling(p/4)*4 - (n1-1)*tildep - barp
			while(addp<0)  addp = addp + 4;
			if(addp>n1-1)  next; 
			if(barp==0)  barP = NULL;
			if(barp>0)  
			{ 
				if(j2==1)  barP = (n2-barp+1):n2
				if(j2==2) {
					if(log2(n2)!=floor(log2(n2)))  next;
					if(log2(n2)>barp)  next; 
					barP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
					if(floor(log2(n2)/2)*2==log2(n2))  barP = c(barP,n2-1)
					if(floor(log2(n2)/2)*2!=log2(n2))  barP = c(barP,2)
					if(length(barP)>barp)  barP = barP[(length(barP)-barp+1):length(barP)]
					if(length(barP)<barp)  barP = c(barP,setdiff(1:n2,barP)[(n2-barp+1):(n2-length(barP))])
				}
				barP = c( sort(setdiff(barP,intersect(barP,tildeP))), sort(intersect(barP,tildeP)) ) 
			}

			groups = NULL  # J: 0 to n1-1; K: 1 to n2.
			groupsadd = NULL
			restK = sort(setdiff(1:n2,tildeP),decreasing=TRUE)
			u = tildeP 
			v = barP
			if(addp>0) if(tildep<n2-2) {
				if(addp>=4) for(k in 1:floor(addp/4)) 
					groupsadd = rbind(groupsadd, c(n1-k*4+(3:0),restK[2],restK[1],restK[1],restK[1]) )
				if(addp - floor(addp/4)*4==3) if(barp>=1) { 
					if(v[1]!=restK[1])  groupsadd = rbind(groupsadd, c(3:0,restK[1],restK[1],restK[1],v[1]) )
					if(v[1]==restK[1])  groupsadd = rbind(groupsadd, c(3:0,restK[2],restK[2],restK[2],v[1]) )
					v = v[-1]
				}
				if(addp - floor(addp/4)*4==3) if(barp==0) { 
					groupsadd = rbind(groupsadd, c(3:1,1,restK[1],restK[1],restK[1],u[1]), 
						c(1,2,3,3,u[c(2,1,1,2)]), c(2,1,2,3,u[c(2,3,3,3)]) )
					u = u[c(-1,-2,-3)]
				}
				if(addp - floor(addp/4)*4==2) if(barp>=2) {
					groupsadd = rbind(groupsadd, c(3:2,0,0,restK[1],restK[1],v[1],v[2]) )
					v = v[c(-1,-2)]
				}
				if(addp - floor(addp/4)*4==2) if(barp==1) { 
					groupsadd = rbind(groupsadd, c(2,1,1,3,restK[1],restK[1],u[1],u[1]), 
						c(2,1,2,3,u[c(1,2,2,2)]), c(0,1,2,3,v[1],u[c(3,3,3)]) )
					u = u[c(-1,-2,-3)]
					v = v[-1]
				}
				if(addp - floor(addp/4)*4==2) if(barp==0) { 
					groupsadd = rbind(groupsadd, c(2,1,1,3,restK[1],restK[1],u[1],u[1]), 
						c(2,1,2,3,u[c(1,2,2,2)]) )
					u = u[c(-1,-2)]
				}
				if(addp - floor(addp/4)*4==1) { 
					theu = u[1] 
					if(length(u)==2) if(length(v)==1) if(u[2]==v[1])  theu = u[2]
					groupsadd = rbind(groupsadd, c(1,1,2,3,restK[2],theu,theu,theu) )
					u = setdiff(u,theu)
				}
			}
			
			inP = sort(intersect(u,v))
			toP = sort(setdiff(u,inP))
			boP = sort(setdiff(v,inP))
			groups = NULL
			if(tildep==n2)  groups = rbind( gpK4(0:(n1-1),barP), gpK4(1:(n1-1),toP) )
			if(tildep==n2-1) if(addp>0)  groups = rbind( gpJ4(1:addp,c(1,tildeP)), gpJ4(setdiff(0:(n1-1),1:addp),tildeP) )
			if(tildep==n2-2) if(addp>0)  groups = rbind( gpJ4(0:addp,c(2,tildeP)), gpJ4(setdiff(0:(n1-1),0:addp),tildeP) )
			if(is.null(groups)) {
				v = c( boP, inP )
				if(length(boP)>0)  u = c( inP, toP )
				if(length(boP)==0) u = c( toP, inP )
				if(length(u)>0) if(length(setdiff(u,v))+length(setdiff(v,u))==0)  u = u[c(length(u),1:(length(u)-1))]
				if(length(v)==0) groups = rbind( gpK4(1:3,u), gpJ4(setdiff(0:(n1-1),0:3),tildeP) )
				if(length(v)>0) if(length(v)<=length(u)) 
					groups = rbind( cbind(0,1,2,3,v,u[1:length(v)],u[1:length(v)],u[1:length(v)]), 
						gpK4(1:3,setdiff(u,u[1:length(v)])), gpJ4(setdiff(0:(n1-1),0:3),tildeP) )
				if(length(v)>length(u)) {
					groups = cbind(0,1,2,3,v[1:length(u)],u,u,u)
					F = matrix(0,((n1/4-1)*tildep),8)
					for(k in 0:((n1/4-2))) for(i in 1:tildep) {
						if(i<tildep-1)  F[tildep*k+i,] = c(4*k+4+0,4*k+4+1,4*k+4+2,4*k+4+3,tildeP[i],tildeP[i],tildeP[i],tildeP[i])
						if(i==tildep-1) F[tildep*k+i,] = c(4*k+4+0,4*k+4+1,4*k+4+0,4*k+4+2,tildeP[i],tildeP[i],tildeP[i+1],tildeP[i+1])
						if(i==tildep  ) F[tildep*k+i,] = c(4*k+4+2,4*k+4+3,4*k+4+1,4*k+4+3,tildeP[i-1],tildeP[i-1],tildeP[i],tildeP[i])
					}
					for(i in 1:((length(v)-length(u))/4))  groups = rbind(groups,
						c(0,0,F[i,1],F[i,2],v[length(u)+4*i-3],v[length(u)+4*i-2],F[i,5],F[i,6]), 
						c(0,0,F[i,3],F[i,4],v[length(u)+4*i-1],v[length(u)+4*i-0],F[i,7],F[i,8]) )
					know = ceiling((length(v)-length(u))/4 /tildep)
					inow= (length(v)-length(u))/4 - (know-1)*tildep
					if(	inow == tildep-1)  groups = rbind( groups, F[(length(v)-length(u))/4+1,] )
					if(	inow <  tildep-1)  groups = rbind( groups, gpJ4(4*know+(0:3),tildeP[(inow+1):tildep]) )
					if(4*know+4 < n1-1)  groups = rbind( groups, gpJ4((4*know+4):(n1-1), tildeP ) )
				}
			}

			if(!is.null(groupsadd)) groups = rbind(groupsadd,groups)
			W = NULL
			for(k in 1:dim(groups)[1]) W = cbind(W, Efull[,groups[k,1:4]*n2+groups[k,5:8]] %*% R4 )
			D = W[,(dim(W)[2]-p+1):dim(W)[2]]*sqrt(85)/16/2 + 1/2
			FinalSep = min(dist(D))
			
			if(FinalSep>BestSep) {
				BestSep = FinalSep
				BestD = D
				Bestn1 = n1
				Bestn2 = n2
				Besttildep = tildep
				Bestbarp = barp
				Bestq = (n1-1)*tildep+barp
			}
		}
	}

	if(BestSep==-1) return(NULL)
	return(BestD)
}



HDOM8 <- function(n,p) {  # Construction of sixteen-level designs
	if(ceiling(p/4)*4>=n) return(NULL)
	BestSep = -1 
	n1list = c(seq(4,n/4,4))
	n1list = n1list[ floor(n / n1list) * n1list == n ]
	n2list = n / n1list
		
	for(indexn in 1:length(n1list)){
		n1 = n1list[indexn]
		n2 = n2list[indexn] 
		if(n2>2) if(floor(n2/4)*4!=n2) next

		A = Hadamard_Matrix(order=n1)
		if(n1==1) A = matrix(1,1,1)
		if(length(A)!=n1^2) next
		for(i in 1:n1) if(A[i,1]==-1) A[i,] = -A[i,]
		
		n2b = n2*2
		while(floor(n2b/16)*16==n2b)  n2b = n2b/2
		if(n2b==8) n2b = 2
		if(n2b==4) n2b = 2
		H = Hadamard_Matrix(order=n2b/2)
		if(n2b>2) if(length(H)==1)  next 
		if(length(H)>1) for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
		while(length(H)!=n2^2)  H = rbind( cbind(H,H), cbind(H,-H) ) 
		Efull = kronecker(A,H)
		
		for(tildep in floor((p-n2+1)/(n1-1)) : ceiling(p/(n1-1)) )  for(j2 in 1:2)  {
			if(tildep < 2) next
			if(tildep > n2) next
			
			j = 2
			if(log2(n2)!=floor(log2(n2))) j = 1
			if(log2(n2)+1>tildep) j = 1
			if(tildep>=n2/4+2) j = 1
			if(j==1)  tildeP = (n2-tildep+1):n2
			if(j==2) {
				tildeP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
				if(floor(log2(n2)/2)*2==log2(n2))  tildeP = c(tildeP,n2-1)
				if(floor(log2(n2)/2)*2!=log2(n2))  tildeP = c(tildeP,2)
				if(length(tildeP)<tildep) tildeP = c(tildeP,setdiff(1:n2,tildeP)[(n2-tildep+1):(n2-length(tildeP))])
				tildeP = sort(tildeP) 
			}
			
			barp = p - (n1-1)*tildep
			if(tildep>=n2-2)  barp = ceiling(max(c(p,(n1-1)*tildep))/4)*4 -(n1-1)*tildep;
			if(barp>min(c(n2-1,(n1-3)*tildep+2)))  barp = min(c(n2-1,(n1-3)*tildep+2));
			if(barp<0)  barp = 0;
			addp = ceiling(p/4)*4 - (n1-1)*tildep - barp
			while(addp<0)  addp = addp + 4;
			if(addp>n1-1)  next; 
			if(barp==0)  barP = NULL;
			if(barp>0)  
			{ 
				if(j2==1)  barP = (n2-barp+1):n2
				if(j2==2) {
					if(log2(n2)!=floor(log2(n2)))  next;
					if(log2(n2)>barp)  next; 
					barP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
					if(floor(log2(n2)/2)*2==log2(n2))  barP = c(barP,n2-1)
					if(floor(log2(n2)/2)*2!=log2(n2))  barP = c(barP,2)
					if(length(barP)>barp)  barP = barP[(length(barP)-barp+1):length(barP)]
					if(length(barP)<barp)  barP = c(barP,setdiff(1:n2,barP)[(n2-barp+1):(n2-length(barP))])
				}
				barP = c( sort(setdiff(barP,intersect(barP,tildeP))), sort(intersect(barP,tildeP)) ) 
			}

			groups = NULL  # J: 0 to n1-1; K: 1 to n2.
			groupsadd = NULL
			restK = sort(setdiff(1:n2,tildeP),decreasing=TRUE)
			u = tildeP 
			v = barP
			if(addp>0) if(tildep<n2-2) {
				if(addp>=4) for(k in 1:floor(addp/4)) 
					groupsadd = rbind(groupsadd, c(n1-k*4+(3:0),restK[2],restK[1],restK[1],restK[1]) )
				if(addp - floor(addp/4)*4==3) if(barp>=1) { 
					if(v[1]!=restK[1])  groupsadd = rbind(groupsadd, c(3:0,restK[1],restK[1],restK[1],v[1]) )
					if(v[1]==restK[1])  groupsadd = rbind(groupsadd, c(3:0,restK[2],restK[2],restK[2],v[1]) )
					v = v[-1]
				}
				if(addp - floor(addp/4)*4==3) if(barp==0) { 
					groupsadd = rbind(groupsadd, c(3:1,1,restK[1],restK[1],restK[1],u[1]), 
						c(1,2,3,3,u[c(2,1,1,2)]), c(2,1,2,3,u[c(2,3,3,3)]) )
					u = u[c(-1,-2,-3)]
				}
				if(addp - floor(addp/4)*4==2) if(barp>=2) {
					groupsadd = rbind(groupsadd, c(3:2,0,0,restK[1],restK[1],v[1],v[2]) )
					v = v[c(-1,-2)]
				}
				if(addp - floor(addp/4)*4==2) if(barp==1) { 
					groupsadd = rbind(groupsadd, c(2,1,1,3,restK[1],restK[1],u[1],u[1]), 
						c(2,1,2,3,u[c(1,2,2,2)]), c(0,1,2,3,v[1],u[c(3,3,3)]) )
					u = u[c(-1,-2,-3)]
					v = v[-1]
				}
				if(addp - floor(addp/4)*4==2) if(barp==0) { 
					groupsadd = rbind(groupsadd, c(2,1,1,3,restK[1],restK[1],u[1],u[1]), 
						c(2,1,2,3,u[c(1,2,2,2)]) )
					u = u[c(-1,-2)]
				}
				if(addp - floor(addp/4)*4==1) { 
					theu = u[1] 
					if(length(u)==2) if(length(v)==1) if(u[2]==v[1])  theu = u[2]
					groupsadd = rbind(groupsadd, c(1,1,2,3,restK[2],theu,theu,theu) )
					u = setdiff(u,theu)
				}
			}
			
			inP = sort(intersect(u,v))
			toP = sort(setdiff(u,inP))
			boP = sort(setdiff(v,inP))
			groups = NULL
			if(tildep==n2)  groups = rbind( gpK4(0:(n1-1),barP), gpK4(1:(n1-1),toP) )
			if(tildep==n2-1) if(addp>0)  groups = rbind( gpJ4(1:addp,c(1,tildeP)), gpJ4(setdiff(0:(n1-1),1:addp),tildeP) )
			if(tildep==n2-2) if(addp>0)  groups = rbind( gpJ4(0:addp,c(2,tildeP)), gpJ4(setdiff(0:(n1-1),0:addp),tildeP) )
			if(is.null(groups)) {
				v = c( boP, inP )
				if(length(boP)>0)  u = c( inP, toP )
				if(length(boP)==0) u = c( toP, inP )
				if(length(u)>0) if(length(setdiff(u,v))+length(setdiff(v,u))==0)  u = u[c(length(u),1:(length(u)-1))]
				if(length(v)==0) groups = rbind( gpK4(1:3,u), gpJ4(setdiff(0:(n1-1),0:3),tildeP) )
				if(length(v)>0) if(length(v)<=length(u)) 
					groups = rbind( cbind(0,1,2,3,v,u[1:length(v)],u[1:length(v)],u[1:length(v)]), 
						gpK4(1:3,setdiff(u,u[1:length(v)])), gpJ4(setdiff(0:(n1-1),0:3),tildeP) )
				if(length(v)>length(u)) {
					groups = cbind(0,1,2,3,v[1:length(u)],u,u,u)
					F = matrix(0,((n1/4-1)*tildep),8)
					for(k in 0:((n1/4-2))) for(i in 1:tildep) {
						if(i<tildep-1)  F[tildep*k+i,] = c(4*k+4+0,4*k+4+1,4*k+4+2,4*k+4+3,tildeP[i],tildeP[i],tildeP[i],tildeP[i])
						if(i==tildep-1) F[tildep*k+i,] = c(4*k+4+0,4*k+4+1,4*k+4+0,4*k+4+2,tildeP[i],tildeP[i],tildeP[i+1],tildeP[i+1])
						if(i==tildep  ) F[tildep*k+i,] = c(4*k+4+2,4*k+4+3,4*k+4+1,4*k+4+3,tildeP[i-1],tildeP[i-1],tildeP[i],tildeP[i])
					}
					for(i in 1:((length(v)-length(u))/4))  groups = rbind(groups,
						c(0,0,F[i,1],F[i,2],v[length(u)+4*i-3],v[length(u)+4*i-2],F[i,5],F[i,6]), 
						c(0,0,F[i,3],F[i,4],v[length(u)+4*i-1],v[length(u)+4*i-0],F[i,7],F[i,8]) )
					know = ceiling((length(v)-length(u))/4 /tildep)
					inow= (length(v)-length(u))/4 - (know-1)*tildep
					if(	inow == tildep-1)  groups = rbind( groups, F[(length(v)-length(u))/4+1,] )
					if(	inow <  tildep-1)  groups = rbind( groups, gpJ4(4*know+(0:3),tildeP[(inow+1):tildep]) )
					if(4*know+4 < n1-1)  groups = rbind( groups, gpJ4((4*know+4):(n1-1), tildeP ) )
				}
			}

			if(!is.null(groupsadd)) groups = rbind(groupsadd,groups)
			W = NULL
			for(k in 1:dim(groups)[1]) W = cbind(W, Efull[,groups[k,1:4]*n2+groups[k,5:8]] %*% R3 )
			D = W[,(dim(W)[2]-p+1):dim(W)[2]]*sqrt(21)/8/2 + 1/2
			FinalSep = min(dist(D))
			
			if(FinalSep>BestSep) {
				BestSep = FinalSep
				BestD = D
				Bestn1 = n1
				Bestn2 = n2
				Besttildep = tildep
				Bestbarp = barp
				Bestq = (n1-1)*tildep+barp
			}
		}
	}
	
	n1 = 2 
	n2 = n / n1 
	if(n2>2) if(floor(n2/4)*4!=n2) if(BestSep==-1) return(NULL)
	if(n2>2) if(floor(n2/4)*4!=n2) return(BestD)


	A = Hadamard_Matrix(order=n1)
	for(i in 1:n1) if(A[i,1]==-1) A[i,] = -A[i,]
		
	n2b = n2*2
	while(floor(n2b/16)*16==n2b)  n2b = n2b/2
	if(n2b==8) n2b = 2
	if(n2b==4) n2b = 2
	H = Hadamard_Matrix(order=n2b/2)
	if(n2b>2) if(length(H)==1) if(BestSep==-1) return(NULL)
	if(n2b>2) if(length(H)==1) return(BestD)
	if(length(H)>1) for(i in 1:dim(H)[1]) if(H[i,1]==-1) H[i,] = -H[i,]
	while(length(H)!=n2^2)  H = rbind( cbind(H,H), cbind(H,-H) ) 
	Efull = kronecker(A,H)
		
	barp = ceiling(p/4)*2
	tildep = ceiling(p/4)*2
	for(j in 1:2)  {
		if(j==2) {
			if(log2(n2)!=floor(log2(n2)))  next;
			if(log2(n2)+1>tildep)  next; 
			tildeP = c(n2 - 2^(1:(log2(n2)-1)) , n2)
			if(floor(log2(n2)/2)*2==log2(n2))  tildeP = c(tildeP,n2-1)
			if(floor(log2(n2)/2)*2!=log2(n2))  tildeP = c(tildeP,2)
			if(length(tildeP)<tildep) tildeP = c(tildeP,setdiff(1:n2,tildeP)[(n2-tildep+1):(n2-length(tildeP))])
			tildeP = sort(tildeP) 
			barP = tildeP
		}
		if(j==1)  { tildeP = (n2-tildep+1):n2;  barP = tildeP; }
		groups = cbind(0,1,0,1,seq(n2-tildep+1,n2-1,2),seq(n2-tildep+1,n2-1,2),seq(n2-tildep+2,n2,2),seq(n2-tildep+2,n2,2))

		W = NULL
		for(k in 1:dim(groups)[1]) W = cbind(W, Efull[,groups[k,1:4]*n2+groups[k,5:8]] %*% R3 )
		D = W[,(dim(W)[2]-p+1):dim(W)[2]]*sqrt(21)/8/2 + 1/2
		FinalSep = min(dist(D))
		
		if(FinalSep>BestSep) {
			BestSep = FinalSep
			BestD = D
			Bestn1 = n1
			Bestn2 = n2
			Besttildep = tildep
			Bestbarp = barp
			Bestq = (n1-1)*tildep+barp
		}
	}

	if(BestSep==-1) return(NULL)
	return(BestD)
}


HDOMdesign <- function(n,p,s) {  # Construction of high-dimensional orthogonal maximin distance designs
	if(s==2) return(HDOM2(n,p))
	if(s==4) return(HDOM4(n,p))
	if(s==8) return(HDOM8(n,p))
	if(s==16) return(HDOM16(n,p))
}
