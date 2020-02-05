from math import log, sqrt, exp
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
#ek = lambda kf, Me: 0.125*dp2*(2.0*kf*Ef(kf,Me)**3-kf*Ef(kf,Me)*Me*Me \
#                               -Me**4*log((Ef(kf,Me)+kf)/Me))
# pk is the kinetic pressure in energy-pressure tensor 
#pk = lambda kf, Me: 1.0/24.0*dp2*(2.0*kf*Ef(kf,Me)**3-5.0*kf*Ef(kf,Me)*Me*Me \
#                               +3.0*Me**4*log((Ef(kf,Me)+kf)/Me))

# Some of the initial value for sigma, omega, rho and Coulumb fields 
global fsig, fomg, frho, fcou

fsig, fomg, frho, fcou = 10., 10., 10., 1.
    
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

    npv, nnv = 0.5*hc3*n*(1.0-a), 0.5*hc3*n*(1.0+a)
    nv = npv+nnv
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me = M - gs*X[0] 
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-f0*(ns(kfp,Me)+ns(kfn,Me))
        Y[1] = X[1]-f1*nv
        return Y  
    fsig, fomg = sf.NewtRoot(X,func)
    
    Me = M - gs*fsig  
    nps, nns = ns(kfp,Me), ns(kfn,Me)
    
    L = 0.5*mw2*fomg**2-0.5*ms2*fsig**2
    Uv = mw2*fomg**2 
    
    ekp = 0.75*npv*Ef(kfp,Me)+0.25*Me*nps
    ekn = 0.75*nnv*Ef(kfn,Me)+0.25*Me*nns
    e = ekp + ekn + Uv - L 
  
    
    pkp = 0.25*npv*Ef(kfp,Me)-0.25*Me*nps
    pkn = 0.25*nnv*Ef(kfn,Me)-0.25*Me*nns
    p = pkp + pkn + L 
    
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

    npv, nnv = 0.5*hc3*n*(1.0-a), 0.5*hc3*n*(1.0+a)
    nv = npv+nnv
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me =  M - gs*X[0]
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-f00*(ns(kfp,Me)+ns(kfn,Me))+f02*X[0]**2+f03*X[0]**3
        Y[1] = X[1]-f10*nv+f13*X[1]**3 
        Y[2] = X[2]-f20*(nnv-npv) 
        return Y
    fsig, fomg, frho = sf.NewtRoot(X,func) 
    
    Me = M - gs*fsig  
    nps, nns = ns(kfp,Me), ns(kfn,Me)
    
    L = 0.5*mw2*fomg**2+0.25*c3*fomg**4+0.5*mr2*frho**2 \
       -0.5*ms2*fsig**2-g2*fsig**3/3.0-0.25*g3*fsig**4
    Uv = mw2*fomg**2 + c3*fomg**4 + mr2*frho**2 
    
    ekp = 0.75*npv*Ef(kfp,Me)+0.25*Me*nps
    ekn = 0.75*nnv*Ef(kfn,Me)+0.25*Me*nns
    e = ekp + ekn + Uv - L 
    
    pkp = 0.25*npv*Ef(kfp,Me)-0.25*Me*nps
    pkn = 0.25*nnv*Ef(kfp,Me)-0.25*Me*nns
    p = pkp + pkn + L 
    
    Ea = 0.0
    if(nv!=0): Ea = e/nv-M
    return Ea, p/hc3, 1.0-gs*fsig/M #, gw*fomg+gr*frho, gw*fomg-gr*frho


def tm1_L40_nmeos(n,a):
    ''' The nuclear matter EOS of the modified TM1 model with L = 40.
        Input n, a as baryon number density and assymetry parameter.
        Output:
            E and P as the energy per nucleon and pressure &
            Us, Uv, the scalar and vector potentials as gmeson*fmeson. 
            All units in the function are related with [MeV].'''
    #  The parameters :
    ms, gs, g2, g3 = 511.198, 10.0289, 7.2325*hc, 0.6183
    mw, gw, c3 = 783.0, 12.6139, 71.3075
    mr, gr = 770.0,  6.9857
    lam = 0.3432*(gr*gw)**2
    M = 938.0
    ms2, mw2, mr2 = ms**2, mw**2, mr**2
    
    f00, f02, f03 = gs/ms2, g2/ms2, g3/ms2
    f10, f12, f13 = gw/mw2, lam/mw2, c3/mw2
    f20, f22 = gr/mr2, lam/mr2 
    
    kfp = (1.0-a)**(1.0/3.0)*kf(n)
    kfn = (1.0+a)**(1.0/3.0)*kf(n)
    
    global fsig, fomg, frho
    X = [fsig, fomg, frho]

    npv, nnv = 0.5*hc3*n*(1.0-a), 0.5*hc3*n*(1.0+a)
    nv = npv+nnv
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me =  M - gs*X[0]
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-f00*(ns(kfp,Me)+ns(kfn,Me))+f02*X[0]**2+f03*X[0]**3
        Y[1] = X[1]-f10*nv+f12*X[2]**2*X[1]+f13*X[1]**3 
        Y[2] = X[2]-f20*(nnv-npv)+f22*X[1]**2*X[2]
        return Y
    fsig, fomg, frho = sf.NewtRoot(X,func) 
    
    Me = M - gs*fsig  
    nps, nns = ns(kfp,Me), ns(kfn,Me)
    
    L = 0.5*mw2*fomg**2+0.25*c3*fomg**4+0.5*mr2*frho**2 \
       -0.5*ms2*fsig**2-g2*fsig**3/3.0-0.25*g3*fsig**4 \
       +0.5*lam*(frho*fomg)**2
    Uv = mw2*fomg**2 + c3*fomg**4 + mr**2*frho**2 + 2.*lam*frho**2*fomg**2
    
    ekp = 0.75*npv*Ef(kfp,Me)+0.25*Me*nps
    ekn = 0.75*nnv*Ef(kfn,Me)+0.25*Me*nns
    e = ekp + ekn + Uv - L 
    
    pkp = 0.25*npv*Ef(kfp,Me)-0.25*Me*nps
    pkn = 0.25*nnv*Ef(kfp,Me)-0.25*Me*nns
    p = pkp + pkn + L 
    
    Ea = 0.0
    if(nv!=0): Ea = e/nv-M
    return Ea, p/hc3, 1.0-gs*fsig/M #, gw*fomg+gr*frho, gw*fomg-gr*frho


def iufsu_nmeos(n,a):
    ''' The nuclear matter EOS of IUFSU.
        Input n, a as baryon number density and assymetry parameter.
        Output:
            E and P as the energy per nucleon and pressure &
            Us, Uv, the scalar and vector potentials as gmeson*fmeson. 
            All units in the function are related with [MeV].'''
    #  The parameters :
    ms, gs, g2, g3 = 491.5, 9.9713, 8.4929*hc, 0.4877
    mw, gw, c3 = 782.5, 13.0321, 144.2195
    mr, gr = 763.0,  6.79495
    lam = 0.368*(gr*gw)**2
    M = 939.0
    ms2, mw2, mr2 = ms**2, mw**2, mr**2
    
    f00, f02, f03 = gs/ms2, g2/ms2, g3/ms2
    f10, f12, f13 = gw/mw2, lam/mw2, c3/mw2
    f20, f22 = gr/mr2, lam/mr2 
    
    kfp = (1.0-a)**(1.0/3.0)*kf(n)
    kfn = (1.0+a)**(1.0/3.0)*kf(n)
    
    global fsig, fomg, frho
    X = [fsig, fomg, frho]

    npv, nnv = 0.5*hc3*n*(1.0-a), 0.5*hc3*n*(1.0+a)
    nv = npv+nnv
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me =  M - gs*X[0]
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-f00*(ns(kfp,Me)+ns(kfn,Me))+f02*X[0]**2+f03*X[0]**3
        Y[1] = X[1]-f10*nv+f12*X[2]**2*X[1]+f13*X[1]**3 
        Y[2] = X[2]-f20*(nnv-npv)+f22*X[1]**2*X[2]
        return Y
    fsig, fomg, frho = sf.NewtRoot(X,func) 
    
    Me = M - gs*fsig  
    nps, nns = ns(kfp,Me), ns(kfn,Me)
    
    L = 0.5*mw2*fomg**2+0.25*c3*fomg**4+0.5*mr2*frho**2 \
       -0.5*ms2*fsig**2-g2*fsig**3/3.0-0.25*g3*fsig**4 \
       +0.5*lam*(frho*fomg)**2
    Uv = mw2*fomg**2 + c3*fomg**4 + mr2*frho**2 + 2.*lam*frho**2*fomg**2
    
    ekp = 0.75*npv*Ef(kfp,Me)+0.25*Me*nps
    ekn = 0.75*nnv*Ef(kfn,Me)+0.25*Me*nns
    e = ekp + ekn + Uv - L 
    
    pkp = 0.25*npv*Ef(kfp,Me)-0.25*Me*nps
    pkn = 0.25*nnv*Ef(kfp,Me)-0.25*Me*nns
    p = pkp + pkn + L 
    
    Ea = 0.0
    if(nv!=0): Ea = e/nv-M
    return Ea, p/hc3, 1.0-gs*fsig/M #, gw*fomg+gr*frho, gw*fomg-gr*frho


def ddme1_nmeos(n,a):
    ''' The nuclear matter EOS of density-dependent meson exchange model.
        Input n, a as baryon number density and assymetry parameter.
        Output:
            E and P as the energy per nucleon and pressure &
            Us, Uv, the scalar and vector potentials as gmeson*fmeson. 
            All units in the function are related with [MeV].'''
    #  The parameters :
    ms, gs, ps = 550.0, 10.72854,(1.365469, 0.226061, 0.409704, 0.901995)
    mw, gw, pw = 783.0, 13.29015,(1.402488, 0.172577, 0.344293, 0.983955)
    mr, gr, ar = 763.0,  3.66098, 0.515 
    
    n0 = 0.153 # the saturation density
    x = n/n0 
    M = 939.0
    
    ms2, mw2, mr2 = ms**2, mw**2, mr**2
    #  The density dependent coupling constant 
    def G(x,a,b,c,d):
        fx = (1.0+b*(x+d)**2)/(1.0+c*(x+d)**2)
        return a*fx 
    n0v = n0*hc3
    #  Derivative of the coupling constant  
    def R(x,a,b,c,d):
        re = 2.0*a*(b-c)*(x+d)/n0v
        return re/(1.0+c*(x+d)**2)**2
            
    fs, fw, fr = G(x,*ps)*gs/ms2, G(x,*pw)*gw/mw2, exp(-ar*(x-1))*gr/mr2  
    
    kfp = (1.0-a)**(1.0/3.0)*kf(n)
    kfn = (1.0+a)**(1.0/3.0)*kf(n)
    
    global fsig, fomg, frho
    X = [fsig, fomg, frho]

    npv, nnv = 0.5*hc3*n*(1.0-a), 0.5*hc3*n*(1.0+a)
    nv = npv+nnv
    n3 = nnv-npv
    # Note 'ns' and 'nv' are in unit [MeV^3]
    def func(X):
        Me =  M - gs*X[0]
        Y = [1.0 for i in range(len(X))]
        Y[0] = X[0]-fs*(ns(kfp,Me)+ns(kfn,Me)) 
        Y[1] = X[1]-fw*nv 
        Y[2] = X[2]-fr*n3
        return Y
    fsig, fomg, frho = sf.NewtRoot(X,func) 
    
    Me = M - gs*fsig*G(x,*ps)  
    nps, nns = ns(kfp,Me), ns(kfn,Me)
    
    # The rearranged terms 
    nR =nv*(-gs*(nps+nns)*fsig*R(x,*ps)+gw*nv*fomg*R(x,*pw)-gr*ar*n3*frho*fr/n0v)
    L = nR - 0.5*ms2*fsig**2+0.5*mw2*fomg**2+0.5*mr2*frho**2
    Uv = mw2*fomg**2+mr2*frho**2
    
    ekp = 0.75*npv*Ef(kfp,Me)+0.25*Me*nps
    ekn = 0.75*nnv*Ef(kfn,Me)+0.25*Me*nns
    e = ekp + ekn + nR + Uv - L 
    
    pkp = 0.25*npv*Ef(kfp,Me)-0.25*Me*nps
    pkn = 0.25*nnv*Ef(kfp,Me)-0.25*Me*nns
    p = pkp + pkn + L 
    
    Ea = 0.0
    if(nv!=0): Ea = e/nv-M
    return Ea, p/hc3, 1.0-gs*fsig/M #, gw*fomg+gr*frho, gw*fomg-gr*frho


def nm__main__(param,alf=0.0,ni=0.0,nf=0.5,dn=0.001):
    ''' Writting for the EOS in given region [ni, nf].
        param is the name of inside parameters, namely:
            walecka model: walecka, w 
                      tm1: tm1, tm. '''
    
    if(param in ('walecka', 'w')):
        nmeos = walecka_nmeos
    elif(param in ('t','tm1', 'tm')):
        nmeos = tm1_nmeos
    elif(param in ('t40','tm1-40', 'L40')):
        nmeos = tm1_L40_nmeos
    elif(param in ('i','iufsu', 'iu-fsu','iu')):
        nmeos = iufsu_nmeos
    elif(param in ('d','d1','ddme','dd1','dd','ddme1')):
        nmeos = ddme1_nmeos
    
    
    n = ni 
    with open(param+'eos.d','w+') as f:
        f.write('   {0:11} {1:12} {2:12} {3:12} \n'.format('n', 'E/A', 'Pres','M*/M'))
        while(n<=nf+dn):
            enea, pres, reff = nmeos(n,alf)
            f.write('{0:7.3f} {1:12.4f} {2:12.4f} {3:12.4f} \n'  \
                    .format(n,enea, pres, reff))
            n += dn
