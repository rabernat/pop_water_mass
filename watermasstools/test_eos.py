import numpy as np
from seawater import eos80
import eos

def test_eos():
    N = 8
    
    shape = (2,4)
    
    
    S = np.linspace(28,40,N)
    T = np.linspace(0,25,N)
    S.shape = shape
    T.shape = shape
    p = 0.
        
    out = eos.state_mod.state(p,T,S)
    jmd95dens0 = out[0]
    alpha = -out[1]/jmd95dens0
    beta = out[2]/jmd95dens0

    print 'Potential Density (p=0dbar)'
    print eos80.dens0(S,T)
    print jmd95dens0
    np.testing.assert_allclose(eos80.dens0(S,T),jmd95dens0, rtol=1e-5)
    
    print 'Alpha'
    print eos80.alpha(S,T,0.)
    print alpha
    np.testing.assert_allclose(eos80.alpha(S,T,0.),alpha, rtol=1e-1)
    
    print 'Beta'
    print eos80.beta(S,T,0.)
    print beta
    np.testing.assert_allclose(eos80.beta(S,T,0.),beta, rtol=1e-3)
    
def main():
    test_eos()

if __name__ == "__main__":
    main()