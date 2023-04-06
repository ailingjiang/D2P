# -*- coding: utf-8 -*-
"""
This python script is written follow by the logic provided in Schwanghart et al.(2017)
with modification tailered for depressions integrated DEM

@author: ailing
"""
import numpy as np
from cvxopt import spmatrix,matrix, sparse,solvers
from Stream_Ordering import streamOrdering
import time
from math import sqrt
import matplotlib.pyplot as plt

def ismember(a,b,concat=True):
    # a: criteria or query
    # b: pool
    Lia = np.isin(a,b)
    sorter = np.argsort(b)
    Locb = []
    r = np.searchsorted(b,a[Lia],side='right',sorter=sorter)
    l = np.searchsorted(b,a[Lia],side='left',sorter=sorter)
    for b, e in zip(l, r):
        Locb.append(sorter[b:e])
    if (concat) & (len(Locb)>0):
        Locb = np.concatenate(Locb)
    return Lia,Locb


def CRS_Carving(dem_refined,depTopRefined, LPC_strseg,LPC_flrdir,tau,K,Kdep, minStrLen,minslope, maxslope):
    ny, nx =  dem_refined.shape
    geotransform = dem_refined.geotransform
    cell_size = np.round(geotransform[1], decimals=2)
    xidx, yidx = np.meshgrid(np.linspace(0, nx-1, nx), np.linspace(0, ny-1, ny))
    dem_RivSmooth = dem_refined.copy()
    S = streamOrdering(LPC_strseg, LPC_flrdir,cell_size)
    
    ##------------------ Calculate minimum elevation on the banks---------------
    minNeiElvDict = {}
    ku = {0:[0,-1],1:[-1,0],2:[0,1],3:[1,0]}
    NoNeiElv = 9999
    for key in S.streamNodes.keys():
        neix = np.zeros((S.streamNodes[key].x.shape[0],4),dtype=int)
        neiy = np.zeros((S.streamNodes[key].x.shape[0],4),dtype=int)
        neiDEM = np.ones((S.streamNodes[key].x.shape[0],4))*NoNeiElv
        for neidir in range(0,4):
            BoolNotStream = np.zeros((S.streamNodes[key].x.shape[0],), dtype=bool)
            BoolNotNeiDep = np.zeros((S.streamNodes[key].x.shape[0],), dtype=bool)
            neix[:,neidir] = S.streamNodes[key].x+ku[neidir][0]
            neiy[:,neidir] = S.streamNodes[key].y+ku[neidir][1]
            BoolNotBound = (neix[:,neidir] >= 0) * (neix[:,neidir] < nx) * (neiy[:,neidir] >= 0) * (neiy[:,neidir] < ny)
            BoolNotStream[BoolNotBound] = LPC_strseg[neiy[BoolNotBound,neidir],neix[BoolNotBound,neidir]] == LPC_strseg.no_data
           # BoolNotNeiDep[BoolNotStream] = (depTopRefined[neiy[BoolNotStream,neidir],neix[BoolNotStream,neidir]] <= 0) + ((depTopRefined[S.streamNodes[key].y[BoolNotStream],S.streamNodes[key].x[BoolNotStream]]> 0) * (depTopRefined[neiy[BoolNotStream,neidir],neix[BoolNotStream,neidir]] > 0)) 
            BoolNotNeiDep[BoolNotStream] = (depTopRefined[neiy[BoolNotStream,neidir],neix[BoolNotStream,neidir]] <= 0)
            neiDEM[BoolNotNeiDep,neidir] = dem_refined[neiy[BoolNotNeiDep,neidir],neix[BoolNotNeiDep,neidir]]

    
        neiDEM = np.where(neiDEM == dem_refined.no_data, NoNeiElv, neiDEM)
        minNeiElvDict[key] = np.min(neiDEM,axis=1)

    ##------------------ CRS Carving--------------------------------------------------
    NoStiffPoints = []
    ListStrID = np.column_stack((S.topoOrder,S.streamId))
    strMerge = []
    strT = time.time()
    for to in range(max(S.topoOrder),0,-1):
        while len(ListStrID[ListStrID[:,0]==to,1]) > 0:
            sortStrLength = np.sort(np.concatenate([S.length[S.streamId==i] 
                                                    for i in ListStrID[ListStrID[:,0]==to,1]]))
            sortListByLength = [x for _, x in sorted(zip(sortStrLength[::-1], 
                                                         ListStrID[ListStrID[:,0]==to,1]))]
            strID = sortListByLength[0]
            strMerge = [strID]
            strlen = S.length[S.streamId==strID].item()   
            strChildlen = 0
            outlet = S.streamNodes[strID].id[0]
            ## Find the downstream reaches and their leaf streams if the stream length is short
            while strlen  < max(3,minStrLen):
                prestrID = S.streamIdc[S.streamId==strID].item()
                if (prestrID == -1):
                    break
                strMerge.append(prestrID)
                strlen = strlen+ S.length[S.streamId==prestrID].item()
                outlet = S.streamNodes[prestrID].id[0]
                strChildList = S.streamId[(S.streamIdc == prestrID) & (S.streamId != strID)]
                for strChild in strChildList:
                    if S.cumlength[S.streamId==strChild]< max(3,minStrLen):
                          strMerge.append(strChild)
                          strChildlen = strChildlen+S.length[S.streamId==strChild].item()
                          strtemp = []
                          strtemp = strtemp+S.streamId[S.streamIdc==strChild].tolist()
                          i = 0
                          while i < len(strtemp):
                              strChildlen = strChildlen+S.length[S.streamId==strtemp[i]].item()
                              strtemp = strtemp+S.streamId[S.streamIdc==strtemp[i]].tolist()
                              i = i+1
                          strMerge=strMerge+strtemp
                strID = prestrID
                
            # Concat the information related to all notes in the (combined) streams
            x = np.concatenate([S.streamNodes[sid].x for sid in strMerge])
            y = np.concatenate([S.streamNodes[sid].y for sid in strMerge])
            reach = np.concatenate([np.repeat(sid, S.length[S.streamId==sid]) for sid in strMerge])
            outDist = np.concatenate([S.streamNodes[sid].outDist for sid in strMerge])
            minNeiElv = np.concatenate([minNeiElvDict[sid] for sid in strMerge],dtype=np.float64)
            idx = np.concatenate([S.streamNodes[sid].id for sid in strMerge])
            idc = np.concatenate([S.streamNodes[sid].idc for sid in strMerge])            

            # Find the last nodes of the upstream reach
            a = np.array(strMerge)
            _, locb = ismember(a,S.streamIdc)
            b = S.streamId[locb]
            lia,_ = ismember(b,a)
            strUp = b[~lia]
                        
            # Number of the stream nodes in the calculation
            strlen = strlen+strChildlen 
            #sys.exit(0)

            # Formulate matrix for CRS algorithm
            I,loc = ismember(idc, idx)
            #------------------    i-i       i     i+1
            CurvIDs = np.column_stack([idc[loc],idc[I],idx[I]])
            CurvIDs = CurvIDs[CurvIDs[:,1]!=outlet]  # remove the curvature constraint for outlet
            #------------------    i-1       i
            GradIDs = np.column_stack([idc[I],idx[I]])
            
            # Extract original elevation and depression info
            elev1 = dem_RivSmooth[y,x]
            strdep = depTopRefined[y,x]
            #depId = np.unique(strdep[strdep>0])
            
            # For plotting purpose only
            elev_orig = dem_refined[y,x]
            elev_dep = elev_orig.copy()
            elev_dep[strdep==depTopRefined.no_data] = np.nan
            
            # Find the spill id of each depression
            depId = np.unique(strdep[strdep>0])
            depSpillid = []
            elevTemp = 99999
            for did in np.arange(0,depId.size):
                dep_downid = np.min(idx[strdep==depId[did]])
                depSpillid.append(idc[idx==dep_downid][0])
                elevmin = elev1[idx==dep_downid][0]
                while elevTemp > elevmin:
                    dep_downid = idc[idx==depSpillid[-1]]
                    elevTemp = elev1[idx==dep_downid]
                    if elevTemp.size > 0:
                        depSpillid.append(dep_downid[0])
                        elevTemp = elevTemp[0]
                    else:
                        elevTemp = elevmin
            depSpillid = np.array(depSpillid)

            # Remove the curvature constraint for defined points (if any)
            CurvIDs = CurvIDs[np.isin(CurvIDs[:,1],NoStiffPoints,invert=True)]
            nCurvRest = CurvIDs.shape[0]

            # Seperate curvature constraint of depression from the rest
            CurvIdxDep = np.isin(CurvIDs[:,1],np.concatenate([idx[strdep>0],depSpillid]))

            # Remove the gradient constraint for depression except the most downstream point
            NoGradPoints = idx[strdep>0]
            
            GradIDs = GradIDs[np.isin(GradIDs[:,1],NoGradPoints,invert=True)]
            nGradRest = GradIDs.shape[0]
            
            # Prepare for the elevation contraints for river segment lower than banks
            minNEDi = minNeiElv[minNeiElv<9999]
            ElevIDs  = idx[minNeiElv<9999]
            nElevRest = ElevIDs.shape[0]
            
            # Convert id matrix to column idx matrix
            CurvColix = np.empty((nCurvRest,3),dtype=np.int32)
            GradColix = np.empty((nGradRest,3),dtype=np.int32)
            ElevColix = np.empty((nElevRest,1),dtype=np.int32)
            
            seqCol = np.arange(strlen,dtype=np.int32)
            Kvec = np.ones((nCurvRest,))*K
            Kvec[CurvIdxDep] = Kdep

            for p in range(0,3):
                _,loc = ismember(CurvIDs[:,p],idx)
                CurvColix[:,p] = seqCol[loc]
                
            for p in range(0,2):
                _,loc = ismember(GradIDs[:,p],idx)
                GradColix[:,p] = seqCol[loc]
                
            if nElevRest > 0:
                _,loc = ismember(ElevIDs,idx)
                ElevColix = seqCol[loc]

            # Check downstream nodes to see if already processed
            downStr = S.streamIdc[S.streamId==reach[idx==outlet]][0]
            nElevDownRest = 0

            if (downStr != -1) and (downStr not in ListStrID[:,1]) and (~np.isin(outlet,NoGradPoints)):
                xdown = S.streamNodes[downStr].x[-1]
                ydown = S.streamNodes[downStr].y[-1]
                elevdown = dem_RivSmooth[ydown,xdown]+abs(minslope)*cell_size 
                nElevDownRest = 1

            # Construct matrix for the higher to processed downstream outlet elev condition
            if nElevDownRest > 0:           
                loc = np.isin(idx,outlet)
                ElevDownColix = seqCol[loc]
                ElevColix = np.array([],dtype=np.int32)
                minNEDi = np.array([])
                nElevRest = 0
                
            # Construct matrix for the lower to upstream outlet elev condition
            if len(strUp) > 0:
                xUp = np.array([S.streamNodes[sid].x[0] for sid in strUp])
                yUp = np.array([S.streamNodes[sid].y[0] for sid in strUp])
                idcUp = np.array([S.streamNodes[sid].idc[0] for sid in strUp])
                elevUp = dem_RivSmooth[yUp,xUp]
                ElevUpIDs = np.unique(idcUp)
                lia,_ = ismember(ElevUpIDs,NoGradPoints)
                ElevUpIDs = ElevUpIDs[~lia]
                minUpElv = np.array([elevUp[idcUp==i].min() for i in ElevUpIDs])
                minUpElv = minUpElv-abs(minslope)*cell_size                             
                nElevUpRest = ElevUpIDs.shape[0]
                
                ElevUpColix = np.empty((nElevUpRest,1),dtype=np.int32)
                _,loc = ismember(ElevUpIDs,idx)
                ElevUpColix = seqCol[loc]
                lia,locb = ismember(ElevColix,ElevUpColix)
                if lia.any():
                    minNEDi[lia] =  np.minimum(minNEDi[lia],minUpElv[locb])
                    ElevColix = np.concatenate((ElevColix, np.delete(ElevUpColix,locb)))
                    minNEDi = np.concatenate((minNEDi, np.delete(minUpElv,locb)))
                    nElevRest = nElevRest+nElevUpRest-locb.size
                else:
                    ElevColix = np.concatenate((ElevColix, ElevUpColix))
                    minNEDi = np.concatenate((minNEDi, minUpElv))
                    nElevRest = nElevRest+nElevUpRest

                del ElevUpIDs, ElevUpColix, minUpElv,
                
            del CurvIDs, GradIDs,ElevIDs              
            
            # Core CRS algorithm
            if (nCurvRest  > 2) & (nGradRest > 1):              
                ## Construct matrix for curvature restriction
                cDij = outDist[CurvColix[:,1]]-outDist[CurvColix[:,0]]
                cDki = outDist[CurvColix[:,2]]-outDist[CurvColix[:,1]]
                cDkj = outDist[CurvColix[:,2]]-outDist[CurvColix[:,0]]
                Svec = Kvec*cell_size**2*sqrt(strlen/(nCurvRest))
                Asd = spmatrix(np.concatenate((np.multiply(Svec,2/(cDij*cDkj)),
                                               np.multiply(Svec,-2/(cDki*cDij)),
                                               np.multiply(Svec,2/(cDki*cDkj)))),
                               np.tile(np.arange(nCurvRest), 3),
                               np.concatenate((CurvColix[:,0],CurvColix[:,1],CurvColix[:,2])))

                if Asd.size[1] < strlen:
                    Asd = sparse([[Asd],[matrix(np.zeros((nCurvRest,strlen-Asd.size[1])))]])
                
                C = sparse([matrix(np.zeros((strlen, strlen))),Asd]);
                
                ## Construct matrix for gradient restriction
                gDij = outDist[GradColix[:,1]]-outDist[GradColix[:,0]]
                Afd = spmatrix(np.concatenate((1/gDij,-1/gDij)), np.tile(np.arange(nGradRest), 2),
                                 np.concatenate((GradColix[:,0],GradColix[:,1])))
                if Afd.size[1] < strlen:
                    Afd = sparse([[Afd],[matrix(np.zeros((nGradRest,strlen-Afd.size[1])))]])
                    
                ## Construct matrix for elevation restriction (elevation lower than stream bank)
                if minNEDi.size > 0:
                    Anei = spmatrix(1/minNEDi, np.arange(nElevRest),ElevColix)                        
                    if Anei.size[1] < strlen:
                        Anei = sparse([[Anei],[matrix(np.zeros((nElevRest,strlen-Anei.size[1])))]])
                        
                ## Construct matrix for elevation restriction (elevation higher then the inlet of processed downstream reach)   
                if nElevDownRest > 0: 
                    AdownElv = spmatrix([-1/elevdown], np.arange(nElevDownRest),ElevDownColix)  
                    if AdownElv.size[1] < strlen:
                        AdownElv = sparse([[AdownElv],[matrix(np.zeros((nElevDownRest,strlen-AdownElv.size[1])))]])

                if minNEDi.size > 0:
                    if nElevDownRest > 0: 
                        A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2)),
                                    matrix(np.zeros((nElevRest, strlen*2))),
                                    matrix(np.zeros((nElevDownRest, strlen*2)))],
                                    [Afd, matrix(np.zeros((strlen*2, strlen))),Anei,AdownElv]])  
                        e = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,)),np.ones((nElevRest,)),-1*np.ones((nElevDownRest,))),axis=0))
                    else:
                        A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2)),
                                    matrix(np.zeros((nElevRest, strlen*2)))],
                                    [Afd, matrix(np.zeros((strlen*2, strlen))),Anei]])  
                        e = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,)),np.ones((nElevRest,))),axis=0))

                elif nElevDownRest > 0: 
                    A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2)),
                                matrix(np.zeros((nElevDownRest, strlen*2)))],
                                [Afd, matrix(np.zeros((strlen*2, strlen))),AdownElv]])  
                    e = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,)),-1*np.ones((nElevDownRest,))),axis=0))
                else:
                    A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2))],
                                [Afd, matrix(np.zeros((strlen*2, strlen)))]])  
                    e = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,))),axis=0))

                     
                H = 2*(C.T*C)
                H = sparse([matrix(np.zeros((strlen*2,strlen*3))),sparse([[matrix(np.zeros((strlen,strlen*2)))], [H]])]); 
                b = matrix(np.concatenate((elev1, np.zeros((nCurvRest,)))));
                c = -2*C.T*b;
                f = matrix([matrix(tau*np.ones((strlen,))),matrix((1-tau)*np.ones((strlen,))),c])
                Aeq = sparse([[spmatrix(1,range(strlen),range(strlen))],
                                     [-spmatrix(1,range(strlen),range(strlen))],[spmatrix(1,range(strlen),range(strlen))]])
                beq = matrix(elev1.astype('float64'))

                if nElevDownRest>0:
                    res = solvers.qp(P=H, q=f, G=A, h=e, A=Aeq, b=beq, kktsolver='chol', 
                                      options={'show_progress':False,'abstol':1e-6,'kktreg':1e-9})
                else:
                    res = solvers.qp(P=H, q=f, G=A, h=e, A=Aeq, b=beq, kktsolver='ldl2', 
                                  options={'show_progress':False,'abstol':1e-6,'kktreg':1e-9})
                elev2 = np.array(res['x'][2*strlen:])
            
            # Quantile Carving
            elif nGradRest > 1:
                f   = matrix(np.concatenate((tau*np.ones((strlen,)),(1-tau)*np.ones((strlen,)),np.zeros((strlen,)))));
                Aeq = sparse([[spmatrix(1,range(strlen),range(strlen))],
                                     [-spmatrix(1,range(strlen),range(strlen))],[spmatrix(1,range(strlen),range(strlen))]])
                beq = matrix(elev1.astype('float64'))
                
                ## Construct matrix for gradient restriction
                gDij = outDist[GradColix[:,1]]-outDist[GradColix[:,0]]
                Afd = spmatrix(np.concatenate((1/gDij,-1/gDij)), np.tile(np.arange(nGradRest), 2),
                                 np.concatenate((GradColix[:,0],GradColix[:,1])))
                if Afd.size[1] < strlen:
                    Afd = sparse([[Afd],[matrix(np.zeros((nGradRest,strlen-Afd.size[1])))]])
                
                ## Construct matrix for elevation restriction (elevation lower than stream bank)
                if minNEDi.size > 0:
                    Anei = spmatrix(1/minNEDi, np.arange(nElevRest),ElevColix)                        
                    if Anei.size[1] < strlen:
                        Anei = sparse([[Anei],[matrix(np.zeros((nElevRest,strlen-Anei.size[1])))]])
                        
                ## Construct matrix for elevation restriction (elevation higher then the inlet of processed downstream reach)   
                if nElevDownRest > 0: 
                    AdownElv = spmatrix([-1/elevdown], np.arange(nElevDownRest),ElevDownColix)  
                    if AdownElv.size[1] < strlen:
                        AdownElv = sparse([[AdownElv],[matrix(np.zeros((nElevDownRest,strlen-AdownElv.size[1])))]])  

                if minNEDi.size > 0:
                    if nElevDownRest > 0:
                        A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2)),
                                     matrix(np.zeros((nElevRest, strlen*2))),matrix(np.zeros((nElevDownRest, strlen*2)))]
                              ,[Afd, matrix(np.zeros((strlen*2, strlen))),Anei,AdownElv]])
                        b = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,)),np.ones((nElevRest,)),-1*np.ones((nElevDownRest,))),axis=0))               
                    else:                     
                        A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2)),
                                     matrix(np.zeros((nElevRest, strlen*2)))]
                              ,[Afd, matrix(np.zeros((strlen*2, strlen))),Anei]])
                        b = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,)),np.ones((nElevRest,))),axis=0))    
                elif nElevDownRest > 0:
                        A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2)),
                                     matrix(np.zeros((nElevDownRest, strlen*2)))]
                                    ,[Afd, matrix(np.zeros((strlen*2, strlen))),AdownElv]])
                        b = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,)),-1*np.ones((nElevDownRest,))),axis=0))    
                else:
                    A = sparse([[matrix(np.zeros((nGradRest, strlen*2))),spmatrix(-1,range(strlen*2),range(strlen*2))]
                          ,[Afd, matrix(np.zeros((strlen*2, strlen)))]])
                    b = matrix(np.concatenate((-minslope*np.ones((nGradRest,)),np.zeros((strlen*2,))),axis=0))    
                    
                res = solvers.lp(c=f,G=A,h=b,A=Aeq,b=beq,solver='glpk')
                elev2 = np.array(res['x'][2*strlen:])
    
            else:
                elev2 = elev1      
            
            elev2Flat = elev2.flatten()
            dem_RivSmooth[y,x] = elev2Flat
            ListStrID = np.delete(ListStrID,np.flatnonzero(np.isin(ListStrID[:,1],strMerge)),0)
            
            if nElevDownRest > 0:
                AdustNei = minNeiElv<elev2Flat
                xx = x[AdustNei]
                yy = y[AdustNei]
                elev2FlatAdj = elev2Flat[AdustNei]
                for neidir in range(0,4):
                    neix = xx+ku[neidir][0]
                    neiy = yy+ku[neidir][1]
                    neix = neix.astype(int)
                    neiy = neiy.astype(int)
    
                    BoolNotBound = (neix >= 0) * (neix < nx) * (neiy >= 0) * (neiy < ny)
                    BoolNotStream = LPC_strseg[neiy,neix] == LPC_strseg.no_data
                    BoolNotDep =  depTopRefined[neiy,neix]<=0
                    BoolAdjut = BoolNotBound*BoolNotStream*BoolNotDep
                    dem_RivSmooth[neiy[BoolAdjut],neix[BoolAdjut]] = np.maximum(elev2FlatAdj[BoolAdjut], dem_RivSmooth[neiy[BoolAdjut],neix[BoolAdjut]])
           
            plt.figure(1000)
            # #path = np.argsort(outDist)
            for s in strMerge:
                rpath = reach == s
                if s == strMerge[0]:
                    plt.plot(outDist[rpath]*cell_size, elev_orig[rpath],'k',label='Before Hydro-conditioned')
                    plt.plot(outDist[rpath]*cell_size, elev_dep[rpath],color='green', 
                              linestyle='dashed', marker='o',markerfacecolor='g',markersize=4,label='Depression nodes')
                    plt.plot(outDist[rpath]*cell_size, elev2[rpath],'r',label='Adapted CRS')
                else:
                    plt.plot(outDist[rpath]*cell_size, elev_orig[rpath],'k')
                    plt.plot(outDist[rpath]*cell_size, elev_dep[rpath],color='green', 
                              linestyle='dashed', marker='o',markerfacecolor='g',markersize=4)
                    plt.plot(outDist[rpath]*cell_size, elev2[rpath],'r')
                plt.axvline(x=outDist[rpath][-1]*cell_size,color = 'grey',linewidth=1,linestyle='-.')
            plt.axvline(x=outDist.min()*cell_size,color = 'grey',linewidth=1,linestyle='-.')
            plt.title('stream id: ' + ','.join(map(str,strMerge)))
            plt.xlabel('Distance Upstream (m)')
            plt.ylabel('Elevation (m)')
            plt.legend()
            #plt.pause
            plt.pause(0.1)
            plt.show()
            plt.clf()
    
    endT = time.time()
    print('time used: %s' %(endT-strT))  
    return dem_RivSmooth