from Spheral1d import *

class PyTensorQ1d(PythonTensorArtificialViscosity1d):

    def __init__(self, Cl, Cq, WT):
        PythonTensorArtificialViscosity1d.__init__(self, Cl, Cq, WT)
        return

    def computeQPiij(self,
                     xi, vi, rhoi, csi,
                     xj, vj, rhoj, csj,
                     etai, etaj,
                     fCli, fCqi, fClj, fCqj):
        vij = vi - vj
        xij = xi - xj
        
        if vij.dot(xij) < 0.0:  # Compression
            mui = vij.dot(etai) / (etai.magnitude2() + self.epsilon2)
            muj = vij.dot(etaj) / (etaj.magnitude2() + self.epsilon2)
            
            Clij = 0.5 * (fCli + fClj) * self.Cl
            Cqij = 0.5 * (fCqi + fCqj) * self.Cq
            
            ei = -Clij * csi * min(0.0, mui) + Cqij * min(0.0, mui)**2
            ej = -Clij * csj * min(0.0, muj) + Cqij * min(0.0, muj)**2
            
            return (Tensor(ei / rhoi), 
                    Tensor(ej / rhoj),
                    rhoi * ei,
                    rhoj * ej)
        else:
            return (Tensor(0.0), Tensor(0.0), 0.0, 0.0)
    
    def label(self):
        return "PyTensorQ"
