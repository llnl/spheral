from Spheral import *

############################################
# 1D
############################################
class TrackConservation1d:
    def __init__(self,
                 db,
                 hydro,
                 verbose=False):
        
        self.db = db
        self.hydro = hydro
        self.verbose = verbose
        self.conservedQuantities=[['cycle','time','linear momentum']]
        
    def periodicWorkFunction(self,cycle,time,dt):
    
        Ptotx = 0.0

        nodeLists = self.db.nodeLists

        mass = self.db.DEMMass
        velocity = self.db.DEMVelocity

        for nodelisti in range(self.db.numNodeLists):
            for i in range(nodeLists[nodelisti].numInternalNodes):
                Ptotx += mass(nodelisti,i)*velocity(nodelisti,i).x
               

        Ptotx = mpi.allreduce(Ptotx,mpi.SUM)

        self.conservedQuantities.append([cycle,time,Ptotx])

        if self.verbose:
            print((" Total Linear Momentum     : %.15f" % Ptotx))

    def deltaLinearMomentumX(self):
        return abs(self.conservedQuantities[-1][2] - self.conservedQuantities[1][2])

############################################
# 2D
############################################
class TrackConservation2d:
    def __init__(self,
                 db,
                 hydro,
                 verbose=False):
        
        self.db = db
        self.hydro = hydro
        self.verbose = verbose
        self.conservedQuantities=[['cycle','time','linear momentum-x','linear momentum-y','rotational momentum-z']]
        
    def periodicWorkFunction(self,cycle,time,dt):
    
        Ptotx = 0.0
        Ptoty = 0.0
        Rtot  = 0.0
        RtotOmega  = 0.0
        RtotVx  = 0.0
        RtotVy  = 0.0

        nodeLists = self.db.nodeLists

        mass = self.db.DEMMass
        radius = self.db.DEMParticleRadius
        velocity = self.db.DEMVelocity
        position = self.db.DEMPosition
        omega = self.hydro.omega

        for nodelisti in range(self.db.numNodeLists):
            for i in range(nodeLists[nodelisti].numInternalNodes):
                Ptotx += mass(nodelisti,i)*velocity(nodelisti,i).x
                Ptoty += mass(nodelisti,i)*velocity(nodelisti,i).y
                RtotOmega += 0.5 * mass(nodelisti,i) * radius(nodelisti,i)**2 * omega(nodelisti,i)
                RtotVx += mass(nodelisti,i) * (-velocity(nodelisti,i).x * position(nodelisti,i).y)
                RtotVy += mass(nodelisti,i) * ( velocity(nodelisti,i).y * position(nodelisti,i).x)

        Ptotx = mpi.allreduce(Ptotx,mpi.SUM)
        Ptoty = mpi.allreduce(Ptoty,mpi.SUM)
        RtotOmega = mpi.allreduce(RtotOmega,mpi.SUM)
        RtotVx = mpi.allreduce(RtotVx,mpi.SUM)
        RtotVy = mpi.allreduce(RtotVy,mpi.SUM)
        Rtot = RtotOmega+RtotVx+RtotVy

        self.conservedQuantities.append([cycle,time,Ptotx,Ptoty,Rtot])

        if self.verbose:
            print((" Total Linear Momentum X     : %.18f" % Ptotx))
            print((" Total Linear Momentum y     : %.18f" % Ptoty))
            print((" Total Rotational Momentum z : %.18f" % Rtot))

    def deltaLinearMomentumX(self):
        return abs(self.conservedQuantities[-1][2] - self.conservedQuantities[1][2])

    def deltaLinearMomentumY(self):
        return abs(self.conservedQuantities[-1][3] - self.conservedQuantities[1][3])
    
    def deltaRotationalMomentumZ(self):
        return abs(self.conservedQuantities[-1][4] - self.conservedQuantities[1][4])



############################################
# 3D
############################################
class TrackConservation3d:
    def __init__(self,
                 db,
                 hydro,
                 verbose=False):
        
        self.db = db
        self.hydro = hydro
        self.verbose = verbose
        self.conservedQuantities=[['cycle','time','linear momentum-x','linear momentum-y','linear momentum-z','rotational momentum-x','rotational momentum-y','rotational momentum-z']]
        
    def periodicWorkFunction(self,cycle,time,dt):
    
        Ptot = Vector3d.zero
        Rtot = Vector3d.zero

        nodeLists = self.db.nodeLists

        mass = self.db.DEMMass
        radius = self.db.DEMParticleRadius
        velocity = self.db.DEMVelocity
        position = self.db.DEMPosition
        omega = self.hydro.omega

        for nodelisti in range(self.db.numNodeLists):
            for i in range(nodeLists[nodelisti].numInternalNodes):
                Ptot += mass(nodelisti,i) * velocity(nodelisti,i)
                Rtot += mass(nodelisti,i) * (0.4 * radius(nodelisti,i)**2 * omega(nodelisti,i) + position(nodelisti,i).cross(velocity(nodelisti,i)))
                
        Ptot = mpi.allreduce(Ptot,mpi.SUM)
        Rtot = mpi.allreduce(Rtot,mpi.SUM)

        self.conservedQuantities.append([cycle,time,Ptot[0],Ptot[1],Ptot[2],Rtot[0],Rtot[1],Rtot[2]])

        if self.verbose:
            print((" Total Linear Momentum X   : %.15f" % Ptot[0]))
            print((" Total Linear Momentum Y   : %.15f" % Ptot[1]))
            print((" Total Linear Momentum Z   : %.15f" % Ptot[2]))
            print((" Rotational Momentum X     : %.15f" % Rtot[0]))
            print((" Rotational Momentum Y     : %.15f" % Rtot[1]))
            print((" Rotational Momentum Z     : %.15f" % Rtot[2]))

    def deltaLinearMomentumX(self):
        return abs(self.conservedQuantities[-1][2] - self.conservedQuantities[1][2])

    def deltaLinearMomentumY(self):
        return abs(self.conservedQuantities[-1][3] - self.conservedQuantities[1][3])

    def deltaLinearMomentumZ(self):
        return abs(self.conservedQuantities[-1][4] - self.conservedQuantities[1][4])
    
    def deltaRotationalMomentumX(self):
        return abs(self.conservedQuantities[-1][5] - self.conservedQuantities[1][5])

    def deltaRotationalMomentumY(self):
        return abs(self.conservedQuantities[-1][6] - self.conservedQuantities[1][6])
    
    def deltaRotationalMomentumZ(self):
        return abs(self.conservedQuantities[-1][7] - self.conservedQuantities[1][7])
