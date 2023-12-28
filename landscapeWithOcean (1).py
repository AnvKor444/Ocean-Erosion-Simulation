import numpy as np

class LandscapeWithOcean(object):

    def __init__(self,NX,NY):

        # These were all global variables in the matlab code
        self.A = np.zeros((NX,NY))

        self.pool     = np.zeros((NX,NY)) #matrix of pooled areas
        self.drain    = np.zeros((NX,NY)) #matrix of draining points
        self.drainage = np.zeros((NX,NY)) #matrix of drainage points(points connecting drains and pools)
        self.ocean    = np.zeros((NX,NY)) #matrix of ocean points

        dt1 = np.dtype('float,int,int')
        self.ZS = np.zeros(NX*NY,dtype=dt1) #Matrix of Z to be sorted.
        self.VOcean = 0.0
        self.AOcean = 0.0
        self.ZBeachLevel = 0.0


    def ComputeOceanVolumeFromOceanLevelParameter(self,Z,NX,NY,oceanLevelParameter):
        Zmin = np.min(Z)
        Zmax = np.max(Z)
        self.ZBeachLevel = Zmin+oceanLevelParameter*(Zmax-Zmin)
        self.VOcean=0.0
        self.AOcean=0
        for x in range(NX):
            for y in range(NY):
                if (Z[x,y]<=self.ZBeachLevel):
                    self.AOcean += 1;
                    self.VOcean += self.ZBeachLevel-Z[x,y];
        print('Minimum elevation          ',Zmin);
        print('Maximum elevation          ',Zmax);
        print('Beach level                ',self.ZBeachLevel);
        print('Ocean volume               ',self.VOcean);
        print('Percentage of ocean surface',self.AOcean/(NX*NY)*100);
        
    def calculate_collection_area(self,Z,NX,NY):

        flowU = np.zeros((NX,NY)) # there is some flow that goes up UP    from (i,j) to (i,j+1)
        flowD = np.zeros((NX,NY)) # there is some flow that goes up DOWN  from (i,j) to (i,j+1)
        flowL = np.zeros((NX,NY)) # there is some flow that goes up LEFT  from (i,j) to (i,j+1)
        flowR = np.zeros((NX,NY)) # there is some flow that goes up RIGHT from (i,j) to (i,j+1)
        
        self.A = np.zeros((NX,NY))
        
        #loop over all cells, from highest to lowest using sorted array ZS
        for i in range(NX*NY-1,-1,-1):
            x = self.ZS[i][1]
            y = self.ZS[i][2]
            
            xU = np.mod(x-1,NX)       # normally x-1 but observe p.b.c.
            xD = np.mod(x+1,NX)       # normally x+1 but observe p.b.c.
            yL = np.mod(y-1,NY)       # normally y-1 but observe p.b.c.
            yR = np.mod(y+1,NY)       # normally y+1 but observe p.b.c.
        
            #compute fraction of flow to all lower cells
            zsum = 0
            if (Z[x,y]>Z[xU,y]):
                zsum = zsum+Z[x,y]-Z[xU,y] 
            if (Z[x,y]>Z[xD,y]):
                zsum = zsum+Z[x,y]-Z[xD,y] 
            if (Z[x,y]>Z[x,yR]):
                zsum = zsum+Z[x,y]-Z[x,yR] 
            if (Z[x,y]>Z[x,yL]):
                zsum = zsum+Z[x,y]-Z[x,yL] 
            if zsum>0:
                if (Z[x,y]>Z[xU,y]):
                    flowU[x,y] = (Z[x,y]-Z[xU,y])/zsum 
                if (Z[x,y]>Z[xD,y]):
                    flowD[x,y] = (Z[x,y]-Z[xD,y])/zsum 
                if (Z[x,y]>Z[x,yR]):
                    flowR[x,y] = (Z[x,y]-Z[x,yR])/zsum 
                if (Z[x,y]>Z[x,yL]):
                    flowL[x,y] = (Z[x,y]-Z[x,yL])/zsum 
            # the first value is rain on the area (multiply by dx^2 outside)
            # I can only receive water from higher cells, for the flow arrays have
            # been computed.
            if self.drain[x,y]:
                self.A[x,y] = 1 + flowL[x,yR]*self.A[x,yR] + flowR[x,yL]*self.A[x,yL] \
                        + flowU[xD,y]*self.A[xD,y] + flowD[xU,y]*self.A[xU,y]

            elif self.pool[x,y]: # if this cell is part of a pool (boundary), 
                            #add all water that flows into this cell to the drainage point
                self.A[self.drainage == self.pool[x,y]] = \
                        self.A[self.drainage ==self.pool[x,y]] \
                       + flowL[x,yR]*self.A[x,yR] + flowR[x,yL]*self.A[x,yL] \
                       + flowU[xD,y]*self.A[xD,y] + flowD[xU,y]*self.A[xU,y]
            elif self.drainage[x,y]: # if this is drainage point of a pool, 
                self.A[x,y] = 1. + self.A[x,y] + flowL[x,yR]*self.A[x,yR] \
                    + flowR[x,yL]*self.A[x,yL] + flowU[xD,y]*self.A[xD,y] \
                    + flowD[xU,y]*self.A[xU,y]


    def pool_check(self,Z,NX,NY):
        # input  Z(current elevation profile)
        # output ZS, drain, pool, drainage
        # global ZS drain pool drainage;

        ###The matrices below have a border of zeros to make checking the boundary easier                                
        # checked = zeros(NX,NY);   #matrix of checked points
        self.pool     = np.zeros((NX,NY)) #matrix of pooled areas
        self.drain    = np.zeros((NX,NY)) #matrix of draining points
        self.drainage = np.zeros((NX,NY)) #matrix of drainage points(points connecting drains and pools)
        self.ocean = np.zeros((NX,NY)) 
        
        # Try a structured array
        self.ZS = np.zeros(NX*NY,dtype='float,int,int') #Matrix of Z to be sorted.

        #copy all NX*NY cell in 1D array ZS 
        for j in range(NY):
            for i in range(NX):
                self.ZS[i+NX*(j)][0] = Z[i,j]
                self.ZS[i+NX*(j)][1] = i
                self.ZS[i+NX*(j)][2] = j

        #self.ZS = sortrows(self.ZS)               #sort the matrix with ascending elevations
        self.ZS = self.ZS[self.ZS.argsort()]

        #print(self.ZS)

        #Set the main draining point, lowest points of all
        xLow = self.ZS[0][1]
        yLow = self.ZS[0][2]
        self.drain[xLow,yLow] = 1
        self.ocean[xLow,yLow] = 1 # begin with lowest elevation as an ocean
        self.ZBeachLevel = Z[xLow,yLow] # elevation at next lowest ocean point (will be updated in the ocean loop)

        sumVOcean   = 0; # sum of ocean volume elements, is zero to begin
        self.AOcean = 1;    # number of ocean cells

        P = 1  #Pool number to identify pools and associated drainages (pool index)

        #Loop over all cell starting the one above the lowest
        for I in range(1,NX*NY):
            x = self.ZS[I][1]
            y = self.ZS[I][2]
            
            if (sumVOcean < self.VOcean):
                self.ocean[x,y] = 1; # set point as part of ocean
                sumVOcean += ( Z[x,y]-self.ZBeachLevel )*self.AOcean;  #add new ocean layer to all existing ocean points
                self.AOcean     += 1;      # increase ocean area
                self.ZBeachLevel = Z[x,y]; # move the beach up!
            else:
                xU = np.mod(x-1,NX)       # normally x-1 but observe p.b.c.
                xD = np.mod(x+1,NX)       # normally x+1 but observe p.b.c.
                yL = np.mod(y-1,NY)       # normally y-1 but observe p.b.c.
                yR = np.mod(y+1,NY)       # normally y+1 but observe p.b.c.
            
                pL = self.pool[x,yL] #Store values for surrounding pools and current maximum elevation of that pool
                pR = self.pool[x,yR]
                pU = self.pool[xU,y]
                pD = self.pool[xD,y]
            
                #if pL and all(all(self.drainage!=pL)):  ###Establish if pool is present and  
                if pL and (self.drainage!=pL).all():  ###Establish if pool is present and  
                    PL = pL
                else:
                    PL = 0

                #if pR and all(all(self.drainage!=pR)):  ###no drainage point already
                if pR and (self.drainage!=pR).all():  ###Establish if pool is present and  
                    PR = pR
                else:
                    PR = 0

                #if pU and all(all(self.drainage!=pU)):  ###exists for that pool (ie. 
                if pU and (self.drainage!=pU).all():  ###Establish if pool is present and  
                    PU = pU
                else:
                    PU = 0

                #if pD and all(all(self.drainage!=pD)):  ###non-draining pool)
                if pD and (self.drainage!=pD).all():  ###Establish if pool is present and  
                    PD = pD
                else:
                    PD = 0
            
                DL = pL and (self.drainage==pL).any()  ###Establish if pool is present and  
                DR = pR and (self.drainage==pR).any()  ###drainage point does already
                DU = pU and (self.drainage==pU).any()  ###exist for that pool (ie. 
                DD = pD and (self.drainage==pD).any()  ###draining pool)
            
                ###Establish if point contacts any draining areas or drainage points
                dL = self.drain[x,yL] or self.drainage[x,yL] or self.ocean[x,yL]
                dR = self.drain[x,yR] or self.drainage[x,yR] or self.ocean[x,yR] 
                dU = self.drain[xU,y] or self.drainage[xU,y] or self.ocean[xU,y]
                dD = self.drain[xD,y] or self.drainage[xD,y] or self.ocean[xD,y]
            
                ####If connected to a pool that does not yet have a drainage point, the
                ####point must either be a pool or a point draining that pool (ie.
                ####drainage)
                if PL or PR or PU or PD:
                    if dL or dR or dU or dD or DL or DR or DU or DD:
                        self.drainage[x,y] = max([PL,PR,PU,PD])
                    else:
                        self.pool[x,y] = max([PL,PR,PU,PD])

                    if PL:
                        self.pool[self.pool==PL] = max([PL,PR,PU,PD])
                    if PR:
                        self.pool[self.pool==PR] = max([PL,PR,PU,PD])
                    if PU:
                        self.pool[self.pool==PU] = max([PL,PR,PU,PD])
                    if PD:
                        self.pool[self.pool==PD] = max([PL,PR,PU,PD])
                  
                ####if connected to a pool with a drainage point it drains freely
                ####to that pool as the drainage point is at a lower elevation.
                ####Also, if the point is touching a drain or drainage point it
                ####also drains.
                elif pL or pR or pU or pD or dL or dR or dU or dD:
                    self.drain[x,y] = 1
                
            ####If the point doesn't contact a pool, drain, or drainage point,
            ####then it must be a new pool.
                else:
                    self.pool[x,y] = P
                    P = P+1 # increase pool index
            #print(ZBeachLevel,self.VOcean,self.AOcean,np.sum(self.ocean))