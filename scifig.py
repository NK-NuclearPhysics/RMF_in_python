import numpy as np 
import matplotlib.pyplot as plt 
import pandas as pd 
import matplotlib 
import matplotlib.patches as patches 
from matplotlib.ticker import AutoMinorLocator 

#------------------#
# Drawing figures: #
#------------------#   
def FigFrame(ax,ticklen_major=8.0,ticklen_minor=4.0, xloc=5,  yloc=5,  labsize=13.5):	
	''' This function has following features
		1. Thicken the frame
		2. Highlight the ticks on both sides

		Function Form:
			FigFrame(ax,ticklen_major,ticklen_minor, xloc,  yloc,  labsize)  

		Varibles:	
		First you should construct the frame --- ``ax'' at first
		Input parameters:
			ax                the frame axises (must be defined)
			ticklen_major     the length of the major ticks 
			ticklen_minor	  the length of the minor ticks 
			xloc              the locator of x axis
			yloc              the locator of y axis
			labsize           the size of the ticks '''
	matplotlib.rcParams['xtick.direction']='in'
	matplotlib.rcParams['xtick.top']=True
	matplotlib.rcParams['ytick.direction']='in'
	matplotlib.rcParams['ytick.right']=True 

	ax.spines['bottom'].set_linewidth(2.0)
	ax.spines['top'].set_linewidth(2.0)
	ax.spines['left'].set_linewidth(2.0)
	ax.spines['right'].set_linewidth(2.0)

	ax.xaxis.set_minor_locator(AutoMinorLocator(xloc))
	ax.yaxis.set_minor_locator(AutoMinorLocator(yloc))

	ax.tick_params(length=ticklen_major,which='major',width=2.0,labelsize=labsize)
	ax.tick_params(length=ticklen_minor,which='minor',width=2.0,labelsize=labsize)
	return None 

def GetData(cfilepass,*n):
    '''The function return the colomns of the file in a given file pass.

    	Function form:
    		GetData(cfilepass,*n)

    	Variables:	
	        cfilepass --- the filepass
	        n = 1, 2, 3, .... --- numbers of the colomn(s) to read.
	            you can input 1, 2, 3, ..., multicolomns will be returned'''
    try:
        A = pd.read_csv(cfilepass, delim_whitespace=True)
    except:
        print('Error, no such file(s).')
        return False
    A = pd.read_csv(cfilepass, delim_whitespace=True)
    A = np.array(A)
    if(len(n)==1):
        i = n[0]
        return A[:,i-1]
    AA = () # define a tuple
    for i in n:
        AA =AA+(A[:,i-1],)
    return AA

#---------------------#
# Numerical Analyser: #
#---------------------#    
def PolInt(x,xx,yy):
	''' The polynomial interpolation 
			x --- the input value satisfying [x0, x1, x, x2, x3]
			xx --- the input array 
			yy --- the input array '''
	n = len(xx)  
	pp = []
	# The polynomial to start with 
	for i in range(n):   # Copy yy
		pp.append(yy[i])

	n = n -1 		
	j = 1 
	while(j<=n):
		for i in range(n-j):
			h = xx[i+j]-xx[i]
			if(abs(h)==0): h = 1e-6 # the minimal interval 
			pp[i] = ((xx[i+j]-x)*pp[i] + (x-xx[i])*pp[i+1])/h
		j+=1
	return pp[0]		

def LagInt(x,xx,yy):
    ''' The Lagrange's inter- and extra- polation:
        given x, NewtPol evaluates the corresponding y 

        Function form:
        	LagInt(x,xx,yy)

        Variables
	        x --- the input point 
	        xx --- the abscissa values 
	        yy --- the function values.  '''

	# If x is one of the input points, return the value.
    for i in range(len(xx)):
        if(xx[i]==x):
            return yy[i]

    Lx = 1.0
    for xi in xx:
        Lx = Lx*(x-xi) 
    dim = len(xx)
    wj = []
    for xj in xx:
        w = 1.0
        for xi in xx:
            if(xi!=xj):
                w = w*(xj-xi)
        dx = x - xj
        wj.append(1.0/w/dx)
    yx = 0.0
    for i in range(dim):
        yx = yx + yy[i]*wj[i] 
    return Lx*yx

def EvalArry(x,xx,yy,Evalf=LagInt):
    ''' locate the input x, and output the evaluate value at the point x.
        Using disection method with Lagrange interpolation.
		
		Function form:
			EvalArry(x,xx,yy,Evalf)

		Variables:
	       	x -- the point to be evaluated, 
	       	xx -- the input x-abscissa
	       	yy -- the input y values 
	       	Evalf -- the interpolation method, the default is LagInt'''
    h = 0.5*abs(xx[1]-xx[0])
    n = len(xx)-1
    if(x>xx[n]+h or x<xx[0]-h):
        print('interpolation beyond range.')
        return False
    ia, ib = 0, n
    # bracket x in the range [ia,ia+1]
    while(ib-ia>1):
        im = int((ia+ib)/2)
        if(xx[im]<=x): 
            ia = im 
        elif(xx[im]>x):
            ib = im 
    if(ia<=2): 
        xi = xx[0:5]
        yi = yy[0:5]
    elif(ia>=n-1):
        xi = xx[n-4:n]
        yi = yy[n-4:n]
    else:
        xi = xx[ia-2:ia+3]
        yi = yy[ia-2:ia+3]
    return Evalf(x,xi,yi)

def FindMin(xx,yy):
    ''' Find the numerical minimum in the array xx and yy.
        return with the abcissa x and the minimal value y.
        	At firt sight, the minimal value xmin must be trapped into xx.

		Function form:
			FindMin(xx,yy)

		Variables:
            xx --- the abscissa in increasing order
            yy --- the values '''
    gold, eps = 0.618034,1e-6 # The golden section and the accuracy
    cgold = 1.0-gold
    n = len(xx)-1
    
    xa, xm, xb = xx[0],cgold*xx[0]+gold*xx[n],xx[n]
    ya, ym, yb = yy[0], EvalArry(xm,xx,yy), yy[n]
    
    # The down hill method， first we assume ym is the minimal
    while(abs(xa-xb)>eps):
        # downhill from the left side
        xg = cgold*xa+gold*xm
        yg = EvalArry(xg,xx,yy)
        # downhill from the left 
        if(yg<ym):
            xb, yb = xm, ym
            xm, ym = xg, yg
        else:
            xa, ya =  xg, yg 
            xm = cgold*xa+gold*xb
            ym = EvalArry(xm,xx,yy)
    xm = cgold*xa+gold*xb
    return  xm,  EvalArry(xm,xx,yy)

def Deriv1(x,xx,yy,eps=0):
    ''' Five-point evaluation of the 1st derivative
        Function form:
            Deriv1(x,xx,yy):
        Variables:
            xx[0:n-1], yy[0:n-1], the input arrays
            x, the point to obtain the evaluation '''
    h = 0.1*(xx[1]-xx[0])
    if(x<xx[0]-h or x>xx[len(xx)-1]+h):
        print('Variable x in the region beyond evaluation.')
        return False
    if(eps!=0): h = eps 
    h2 = 2*h

    x0, x1, x3, x4 = x - h2, x - h, x+h, x+h2
    dy1 = EvalArry(x3,xx,yy)-EvalArry(x1,xx,yy)
    dy2 = EvalArry(x4,xx,yy)-EvalArry(x0,xx,yy)
    return (8*dy1-dy2)/6.0/h2

def Deriv2(x,xx,yy,eps=0):
    ''' Five-point evaluation of the 2nd derivative
        Function form:
            Deriv2(x,xx,yy):
        Variables:
            xx[0:n-1], yy[0:n-1], the input arrays
             x --- the point to obtain the evaluation 
            eps--- the accuracy, default 1/1000
            
            WARNNING: the calculation on the grids 
                      is not right with NewtPol !!!'''
    h = 0.1*(xx[1]-xx[0])
    if(x<xx[0]-h or x>xx[len(xx)-1]+h):
        print('Variable x in the region beyond evaluation.')
        return False
    
    if(eps!=0): h = eps

    h2 = 2*h
    x0, x1, x3, x4 = x-h2, x-h, x+h, x+h2   
 
    sy1 = EvalArry(x3,xx,yy)+EvalArry(x1,xx,yy)
    sy2 = EvalArry(x4,xx,yy)+EvalArry(x0,xx,yy)
    sy0 = EvalArry(x,xx,yy) 

    return (16.0*sy1 -sy2 - 30*sy0)/12.0/h/h

#-----------------#
# Find the roots: #
#-----------------#
def DisecRoot(xa,xb,Func,eps=1e-6):
    ''' Disection method to search for the 
        roots in given region [xa, xb].
        Function form:
            DisRoot(xa,xb,Func)
        Variables:
            xa, xb --- the searching region for 
                       the root(s), 
            Func --- the input function 
            eps --- the accuracy, default 1e-5'''
    if(xa>xb): # exchange [xa,xb] if xa>xb 
        t = xa 
        xb = xa 
        xa = t 
    xm = 0.5*(xa+xb)
    ya, ym, yb = Func(xa),Func(xm),Func(xb)
    while(xb-xa>=eps):
        if(ya*ym>0 and ym*yb<=0):
            ya = ym
            xa = xm
        elif(ya*ym<=0 and ym*yb>0):
            yb = ym
            xb = xm
        else:
            # the root is not bracket here
            print('the root is not bracket here.')
            return False 
        xm = 0.5*(xa+xb)
        ym = Func(xm)
    return xm 
   

def NewtRoot(X,Func,eps = 1e-6):
    ''' Applying Newton's method to search the root
        of the multidimensional function 
        Function Form:
            NewtRoot(X,Func)
        Variables:
            X --- the multidimensional numpy array
            Func ---- the function user defined
            eps --- the accuracy 
        Note: Func should return with Y in the same 
              dimension as X. '''
    from numpy import abs, sum, eye, max,min
    n = len(X)
    dX = [1.0 for i in range(n)]
    Y= [1.0 for i in range(n)]
    JA = eye(n)
    #JA=array([[1.0 for i in range(n)] for j in range(n)])
    # Note: we can not generate JA by [dX --colomn]
    iterat = 1 
    while(max(abs(Y))>eps):
        Y = Func(X)
        h = 0.1*(min(abs(dX)))
        if(h==0): h = 0.01
        for i in range(n):# The forward derivatives 
            X[i]+= h
            Yh = Func(X)
            X[i]-= h
            for j in range(n):
                JA[j][i]=Yh[j]   

        dX = LUSolv(JA,Y) 
        tot = sum(dX)
        if(tot==1): tot = 1.0+eps 
        for i in range(n):
            dX[i] = h/(-1.0+tot)*dX[i]
            X[i] = X[i] + dX[i]
        iterat += 1 
        if(iterat>=100): return False 
    return  X


#--------------------#
# Matrix operations: #
#--------------------#
def LUSolv(AA,B):
    ''' LU decomposition for the input matrix AA,
        The returning depends on B:
         
    '''
    from numpy import eye
    n = len(AA[0])
    LU =  eye(n)# LU decomposed one 

    for i in range(n):
        LU[i][:] = AA[i][:]
        # The lower part 
        for j in range(i):
            tot = AA[i][j] 
            # For immediate summation 
            for k in range(j):
                tot -= LU[i][k]*LU[k][j]
            if(LU[j][j]==0):LU[j][j] = 1e-6
            LU[i][j]= tot/LU[j][j]
        # The upper part 
        for j in range(i,n):
            LU[i][j] =AA[i][j]
            for k in range(i):
                LU[i][j] -= LU[i][k]*LU[k][j]
    # Outputs depends on the type of B
    tB = type(B)
    if tB is int:
        if(B == 0):# Return the determinant 
            d = 1.0
            for i in range(n):
                d = d*LU[i][i]
            return d
        elif(B == 1): 
            return LU 
        elif(B == -1):  
        # Return RA -- the inversed matrix of AA. 
            RU, RL = eye(n), eye(n)
            # For the upper part of the insersed LU
            for i in range(n-1,-1,-1):
                d = LU[i][i]
                try:
                    RU[i][i] = 1.0/d 
                except:
                    print('Error, singular matrix encountered.')
                    return False
                for j in range(i+1,n):
                    tot = 0.0
                    for k in range(i+1,j+1):
                        tot-=LU[i][k]*RU[k][j] 
                    RU[i][j] = tot/d 
                l = n-1-i
                for j in range(l):
                    tot = -LU[l][j]
                    for k in range(j+1,l):
                        tot -= LU[l][k]*RL[k][j]
                    RL[l][j] = tot
            RA = MatMul(RU,RL)
            return RA
        else:
            return False 
    elif tB is list:
        #B = array(B)
        X = [1.0 for i in range(n)]
        for i in range(n):
            tot = B[i]
            for j in range(i):
                tot -= LU[i][j]*X[j]
            X[i] = tot 
        for i in range(n-1,-1,-1):
            tot = X[i]
            for j in range(n-1-i):
                tot -= LU[i][n-1-j]*X[n-1-j]
             
            if(LU[i][i]==0): LU[i][i] = 1e-6 
            X[i] = tot/LU[i][i]
             
        return X
    
def MatMul(X,Y):
    ''' The matrix multiplication
        Function form:
            MatMul(X,Y)
        Inputs:
            dimensions: 
            (1) dim X, Y = 0, 0
            (2) dim X, Y = 1, 0 or 0, 1 
            (3) dim X, Y = 1, 1
            (4) dim X, Y = 2, 0 or 0, 2
            (5) dim X, Y = 2, 1
            (6) dim X, Y = 2, 2 
            (.) Exceptions will be wrong
        Output:
            (1), (2)  Z = X*Y
            (3)   Z = X dot Y  
            (4)       Z = X*Y
            (5), (6)  Z = X Y '''
    from numpy import array
    dimx, dimy = array(X).ndim, array(Y).ndim
    if(dimx==0 and dimy==0): 
        Z = 1.0*X*Y
    elif(dimx == 1 and dimy == 0):
        Z = [ x*Y for x in X]
    elif(dimx == 0 and dimy == 1):
        Z = [ y*X for y in Y]
    elif(dimx ==1 and dimy == 1):
        n = len(X)
        Z = 0.0 
        for i in range(n):
            Z+= X(i)*Y(i)
    elif(dimx == 2 and dimy == 0):
        Z = [[1.0*x*Y for  x in Ax] for Ax in X]
    elif(dimx == 0 and dimy ==2 ):    
        Z = [[1.0*y*X for  y in Ay] for Ay in Y] 
    elif(dimx == 2 and dimy == 1):
        n = len(Y)
        # Start a new vector 
        Z = [0.0 for i in range(n)]
        for i in range(n):
            for j in range(n):
                Z[i] += X[i][j]*Y[j]
    elif(dimx==2 and dimy==2):
        n = len(Y[0])
        # Start a zero matrix nxn
        Z = [[0.0 for i in range(n)] for j in range(n)]
        for i in range(n):
            for j in range(n):
                for k in range(n):
                    Z[i][j] += X[i][k]*Y[k][j] 
    else:
        return False
    return Z