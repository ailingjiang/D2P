# -*- coding: utf-8 -*-
"""
Created on Tue Mar 30 02:23:37 2021

@author: ailing
"""
import numpy as np
from itertools import compress
## dot dict class is obtained and adapted from the comments in
# https://stackoverflow.com/questions/2352181/how-to-use-a-dot-to-access-members-of-dictionary
# https://gist.github.com/miku/dc6d06ed894bc23dfd5a364b7def5ed8

class dotdict(dict):     
    """dot.notation access to dictionary attributes"""      
    def __getstate__(self):
        state = self.__dict__.copy()
        return state
    def __setstate__(self,state):
        self.__dict__.update(state)

    def __getattr__(*args):         
        val = dict.get(*args) 
        if type(val) is dict:
            return dotdict(val)
        elif val is None:
            raise KeyError(val)
        else:
            return val      
    __setattr__ = dict.__setitem__     
    __delattr__ = dict.__delitem__  


class StreamOrder:
    def __init__(self):
        self.streamId = []
        self.streamIdc = []
        self.topoOrder = []
        self.length = []
        self.cumlength = []
        self.streamNodes = {}
        
    def addFirstNode(self, streamId, id, x, y, idc, outDist):
        self.streamNodes[streamId] = dotdict({'id':[id], 'x':[x], 'y':[y], 'idc':[idc], 'outDist':[outDist]})
        #self.streamNodes[streamId] = {'id':[id], 'x':[x], 'y':[y], 'idc':[idc], 'outDist':[outDist]}
        
    def listOrNpyAppend(self,obj,var):
        if isinstance(obj,list):
            obj.append(var)
        else:
            np.append(obj,var)
    def addNodes(self, streamId, id, x, y, idc, outDist):
        self.listOrNpyAppend(self.streamNodes[streamId].id,id)
        self.listOrNpyAppend(self.streamNodes[streamId].x,x)
        self.listOrNpyAppend(self.streamNodes[streamId].y,y)
        self.listOrNpyAppend(self.streamNodes[streamId].idc,idc)
        self.listOrNpyAppend(self.streamNodes[streamId].outDist,outDist)
        # self.listOrNpyAppend(self.streamNodes[streamId]['id'],id)
        # self.listOrNpyAppend(self.streamNodes[streamId]['x'],x)
        # self.listOrNpyAppend(self.streamNodes[streamId]['y'],y)
        # self.listOrNpyAppend(self.streamNodes[streamId]['idc'],idc)
        # self.listOrNpyAppend(self.streamNodes[streamId]['outDist'],outDist)
                               
    def addStreams(self, streamId, streamIdc, topoOrder):
        self.listOrNpyAppend(self.streamId,streamId)
        self.listOrNpyAppend(self.streamIdc,streamIdc)
        self.listOrNpyAppend(self.topoOrder,topoOrder)
        
    def addLength(self, length, cumlenth):
        self.listOrNpyAppend(self.length,length)
        self.listOrNpyAppend(self.cumlength,cumlenth)
        
        
    def list2npy(self): 
        if isinstance(self.streamId,list):
            self.streamId = np.array(self.streamId)
        if isinstance(self.streamIdc,list):
            self.streamIdc = np.array(self.streamIdc)
        if isinstance(self.topoOrder,list):
             self.topoOrder = np.array(self.topoOrder)
        if isinstance(self.length,list):
            self.length = np.array(self.length)
        if isinstance(self.cumlength,list):
             self.cumlength = np.array(self.cumlength)
                
        for key in list(self.streamNodes.keys()):
            if isinstance(self.streamNodes[key].id,list):
                self.streamNodes[key].id = np.array(self.streamNodes[key].id)
            if isinstance(self.streamNodes[key].x,list):
                self.streamNodes[key].x = np.array(self.streamNodes[key].x)
            if isinstance(self.streamNodes[key].y,list):
                self.streamNodes[key].y = np.array(self.streamNodes[key].y)
            if isinstance(self.streamNodes[key].idc,list):
                self.streamNodes[key].idc = np.array(self.streamNodes[key].idc)
            if isinstance(self.streamNodes[key].outDist,list):
                self.streamNodes[key].outDist = np.array(self.streamNodes[key].outDist)

            # if isinstance(self.streamNodes[key]['id'],list):
            #     self.streamNodes[key]['id'] = np.array(self.streamNodes[key]['id'])
            # if isinstance(self.streamNodes[key]['x'],list):
            #     self.streamNodes[key]['x'] = np.array(self.streamNodes[key]['x'])
            # if isinstance(self.streamNodes[key]['y'],list):
            #     self.streamNodes[key]['y'] = np.array(self.streamNodes[key]['y'])
            # if isinstance(self.streamNodes[key]['idc'],list):
            #     self.streamNodes[key]['idc'] = np.array(self.streamNodes[key]['idc'])
            # if isinstance(self.streamNodes[key]['outDist'],list):
            #     self.streamNodes[key]['outDist'] = np.array(self.streamNodes[key]['outDist'])
            
def streamOrdering(LPC_strseg, LPC_flrdir,cell_size):      
    ny, nx =  LPC_strseg.shape    
    xidx, yidx = np.meshgrid(np.linspace(0, nx-1, nx), np.linspace(0, ny-1, ny))
    outlet_x = np.ma.masked_array(xidx,(LPC_flrdir>=0) | (LPC_strseg==LPC_strseg.no_data) ,fill_value=-1).compressed()
    outlet_y = np.ma.masked_array(yidx,(LPC_flrdir>=0) | (LPC_strseg==LPC_strseg.no_data) ,fill_value=-1).compressed()
    outstream = np.ma.masked_array(LPC_strseg,(LPC_flrdir>=0) | (LPC_strseg==LPC_strseg.no_data) ,fill_value=-1).compressed()
    ## r.watershed direction integer mapping
    ku = {2:[0,-1],4:[-1,0],6:[0,1],8:[1,0]} # python cell offset based on flow direction
    direction = {2:[4,6,8],4: [2,6,8],6:[2,4,8],8:[2,4,6]} # upstream directions
    recRelation = {2:6,4:8,6:2,8:4} #'give and receive' relationship
    
    S = StreamOrder()
    NodeId = -1
    # start with the most downstream reaches
    for outstr in outstream:
        xic = -1; yic = -1; idc = -1; dist = 0;
        NodeId = NodeId+1
        xi = int(outlet_x[outstream==outstr])
        yi = int(outlet_y[outstream==outstr])
        streamTributary =  {outstr:[xi,yi,-1,1,NodeId,idc,dist]}
        # if there is stream not processed in stream triburary
        while len(streamTributary) > 0:       
            stream = list(streamTributary.keys())[0]
            [xi, yi] = streamTributary[stream][:2]       
            S.addStreams(stream,streamTributary[stream][2],streamTributary[stream][3])
            S.addFirstNode(stream,streamTributary[stream][4],xi,yi,
                           streamTributary[stream][5],streamTributary[stream][6])
            idc = streamTributary[stream][4]
            dist = streamTributary[stream][6]
            # if there is upstream cell
            while (xi != xic) | (yi != yic): 
                xic = xi
                yic = yi
                curfdir = abs(LPC_flrdir[yic,xic])
                dist = dist+cell_size
                # search for upstream directions
                for d in direction[curfdir]:
                    tempx = int(xic+ku[d][0])
                    tempy = int(yic+ku[d][1])
                    # if two cells have 'give and receive' relationship and if it is a stream cell
                    if (LPC_strseg[tempy,tempx] != LPC_strseg.no_data) & (recRelation[d] == abs(LPC_flrdir[tempy,tempx])):
                        if LPC_strseg[tempy,tempx] == stream:
                            xi = tempx
                            yi = tempy                                       
                            NodeId = NodeId+1
                            S.addNodes(stream,NodeId,xi,yi,idc,dist) 
                            idcUpdate = NodeId
                        else:
                            NodeId = NodeId+1 
                            streamTributary[LPC_strseg[tempy,tempx]]= [tempx,tempy,
                                                                       stream,streamTributary[stream][3]+1,
                                                                       NodeId,idc,dist]                            
                    # if LPC_strseg[tempy,tempx] == stream:
                    #     xi = tempx
                    #     yi = tempy                                       
                    #     NodeId = NodeId+1
                    #     S.addNodes(stream,NodeId,xi,yi,idc,dist) 
                    #     idcUpdate = NodeId
                    # elif LPC_strseg[tempy,tempx] != LPC_strseg.no_data:
                    #     # if two cells have 'give and receive' relationship
                    #     if recRelation[d] == abs(LPC_flrdir[tempy,tempx]):
                    #         NodeId = NodeId+1 
                    #         streamTributary[LPC_strseg[tempy,tempx]]= [tempx,tempy,
                    #                                                    stream,streamTributary[stream][3]+1,
                    #                                                    NodeId,idc,dist]
                idc = idcUpdate                      
            streamTributary.pop(stream, None)
    
    leafstreams = list(set(S.streamIdc).symmetric_difference(S.streamId))
    lendict = {}
    i = 0
    while i < len(leafstreams):
        if leafstreams[i] != -1:
            leng =  len(S.streamNodes[leafstreams[i]].id)
            booltemp = S.streamIdc==leafstreams[i]
            if sum(booltemp) == 0:
                childlen = 0
            else:
                childstr = list(compress(S.streamId, booltemp))
                childdict = {k: lendict[k] for k in lendict.keys() & childstr}
                childlen = max(int(d['cumlength']) for d in childdict.values())
                
            lendict[leafstreams[i]] = {'length': leng,'cumlength':leng+childlen}
            leafstreams.append(S.streamIdc[S.streamId.index(leafstreams[i])])
        i = i+1
        
    for strId in S.streamId:
        S.addLength(lendict[strId]['length'],lendict[strId]['cumlength'])  
           
        
    S.list2npy()
    return S