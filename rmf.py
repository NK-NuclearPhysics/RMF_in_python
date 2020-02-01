from math import log, sqrt 
import scifig as sf 
import sys 
# definition of the parameters:
hc = 197.327
hc3 = hc**3
pi = 3.141592653589793
dp2 = 1.0/pi/pi


# definition of some functions:
kf = lambda n: (1.5*n/dp2)**(1.0/3.0)*hc

Ef = lambda kf, Me: sqrt(kf**2+Me**2)
# The scalar density of baryon
ns = lambda kf, Me: 0.5*Me*dp2*(kf*Ef(kf,Me)-Me**2*log((kf+Ef(kf,Me))/abs(Me)))
# ek is the kinetic energy density in energy-pressure tensor 
ek = lambda kf, Me: 0.125*dp2*(2.0*kf*Ef(kf,Me)**3-kf*Ef(kf,Me)*Me*Me \
                               -Me**4*log((Ef(kf,Me)+kf)/Me))
# pk is the kinetic pressure in energy-pressure tensor 
pk = lambda kf, Me: 1.0/24.0*dp2*(2.0*kf*Ef(kf,Me)**3-5.0*kf*Ef(kf,Me)*Me*Me \
                               +3.0*Me**4*log((Ef(kf,Me)+kf)/Me))

# Some of the initial value for sigma, omega, rho and Coulumb fields 
global fsig, fomg, frho, fphn

fsig, fomg, frho, fphn = 10., 10., 10., 1.
    
def walecka_nmeos(n,a):
    ''' The nuclear matter EOS of Walecka model.
        Input n, a as baryon number density and assymetry parameter.
        Output:
            E and P as the energy per nucleon and pressure &
            Us, Uv, the scalar and vector potentials as gmeson*fmeson. 
            All units in the function are related with [MeV].'''
    # The parameters : 
    ms, gs = 550.0,  9.5726 
    mw, gw = 783.0, 11.6711
    M = 939.0
    ms2, mw2 = ms**2, mw**2
    f0 = gs/ms2 
    f1 = gw/mw2
        
    kfp = (1.0-a)**(1.0/3.0)*kf(n)
    kfn = (1.0+a)**(1.0/3.0)*kf(n)

    global fsig, fomg
    X = [fsig, fomg]

    nv = hc3*n 
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me = M - gs*X[0] 
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-f0*(ns(kfp,Me)+ns(kfn,Me))
        Y[1] = X[1]-f1*(kfp**3+kfn**3)/3.0*dp2
        return Y  
    fsig, fomg = sf.NewtRoot(X,func)
    
    Me = M - gs*fsig  
    vs = 0.5*ms2*fsig**2
    vw = 0.5*mw2*fomg**2
    e = ek(kfp,Me)+ek(kfn,Me)+vs+vw
    p = pk(kfp,Me)+pk(kfn,Me)-vs+vw
    
    Ea = 0.0
    if(nv!=0): Ea = e/nv-M
 
    return Ea, p/hc3, 1.0-gs*fsig/M  #, gw*fomg 

def tm1_nmeos(n,a):
    ''' The nuclear matter EOS of TM1 model.
        Input n, a as baryon number density and assymetry parameter.
        Output:
            E and P as the energy per nucleon and pressure &
            Us, Uv, the scalar and vector potentials as gmeson*fmeson. 
            All units in the function are related with [MeV].'''
    #  The parameters :
    ms, gs, g2, g3 = 511.198, 10.0289, 7.2325*hc, 0.6183
    mw, gw, c3 = 783.0, 12.6139, 71.3075
    mr, gr = 770.0,  4.6322 
    M = 938.0
    ms2, mw2, mr2 = ms**2, mw**2, mr**2
    
    f00, f02, f03 = gs/ms2, g2/ms2, g3/ms2
    f10, f13 = gw/mw2, c3/mw2
    f20 = gr/mr2
    
    kfp = (1.0-a)**(1.0/3.0)*kf(n)
    kfn = (1.0+a)**(1.0/3.0)*kf(n)
    
    global fsig, fomg, frho
    X = [fsig, fomg, frho]
    
    nv = hc3*n 
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me =  M - gs*X[0]
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-f00*(ns(kfp,Me)+ns(kfn,Me))+f02*X[0]**2+f03*X[0]**3
        Y[1] = X[1]-f10*nv+f13*X[1]**3 
        Y[2] = X[2]-f20*a*nv 
        return Y
    fsig, fomg, frho = sf.NewtRoot(X,func) 
 
    Me = M - gs*fsig
    vs = 0.5*ms2*fsig**2+g2*fsig**3/3.0+0.25*g3*fsig**4
    vw = 0.5*mw2*fomg**2+0.25*c3*fomg**4
    vr = 0.5*mr2*frho**2
    
    e = ek(kfp,Me)+ek(kfn,Me)+vs+vw+vr+0.5*c3*fomg**4
    p = pk(kfp,Me)+pk(kfn,Me)-vs+vw+vr
    
    Ea = 0.0
    if(nv!=0): Ea = e/nv-M
    return Ea, p/hc3, 1.0-gs*fsig/M #, gw*fomg+gr*frho, gw*fomg-gr*frho



def nm__main__(param,alf=0.0,ni=0.0,nf=0.5,dn=0.001):
    ''' Writting for the EOS in given region [ni, nf].
        param is the name of inside parameters, namely:
            walecka model: walecka, w 
                      tm1: tm1, tm. '''
    
    if(param in ('walecka', 'w', 'W', 'Walecka')):
        nmeos = walecka_nmeos
    elif(param in ('t','tm1', 'tm', 'TM','TM1')):
        nmeos = tm1_nmeos
    
    n = float(ni)
    with open(param+'eos.d','w+') as f:
        f.write('   {0:11} {1:12} {2:12} {3:12} \n'.format('n', 'E/A', 'Pres','M*/M'))
        while(n<=float(nf)+float(dn)):
            enea, pres, reff = nmeos(n,float(alf))
            f.write('{0:7.3f} {1:12.4f} {2:12.4f} {3:12.4f} \n'  \
                    .format(n,enea, pres, reff))
            n += dn

''' Suggestions:
		Input parameters:
            Walecka:  w, walecka, W, Walecka
            TM1: TM, TM1, tm1, tm1
'''
            
cname = sys.argv[1:]
nm__main__(*cname)    