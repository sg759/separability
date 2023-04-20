#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun May  1 20:49:59 2022

"""

import sys
import numpy as np
# import pdb

class read_lmp_file(object):
    
    def __init__(self,filename, filetype = 'lammps', using_mex = 'no', tstep = -1):
        self.inputFileType = filetype
        self.filename = filename
        self.using_mex = using_mex 
        # self.nratoms
        
        # mex_specified = 0
        
        ## default values 
        self.atom_style = 'full'
        self.nrtypes = 0
        self.periodicity = [{'pp'}, {'pp'} ,{'pp'}]
        self.velocity = np.empty(shape=(0,3))
        self.molecule = np.empty(shape=(0,1)) 
        self.charge = np.empty(shape=(0,1))
        self.image_flag = np.empty(shape=(0,3)) 
        self.added = np.empty(shape=(0,0))
        self.type = np.empty(shape=(0,1)) 
        self.nratoms = 0
        self.timestep = -1
        self.fields = []
        self.extra_fields = []
        self.write_extras = 0
        self.extras = []
        self.errorflag = 0
        

    def read_lmp_file(self):
        
        # 'filetype'  :   string 'lammps' or 'dump'  [lammps]
        # timestep'  :   integer 
        
        # print ('USING lmp_file_cleaned version \n')
        
        if self.using_mex == 'no':
            self.fid = self.open_file()
        tstep = -1
        
        if self.inputFileType == 'lammps':
            self = self.load_data(tstep)
        # self.nrtypes = max(self.type)
        if self.inputFileType == 'dump':
            self = self.next_step(tstep)

        
        return self
    
    
    def next_step(self,tstep = -1):
        ''' read a single step in a dump file'''
        
        self.fields = []
        
        if not self.filename:
            sys.stderr.write( "lmp_class.py:next_step - Error: Oject filename is not defined\n" )
            sys.exit(-1)
        elif self.fid == -1:
            self.file_open()
      

      # --------------------- default values
      
        correct_timestep = -1
        if tstep != -1: 
            correct_timestep = 1
       
        present = []
        ttest = self.timestep;            # check for first call... memory allocation
        nratoms_old = self.nratoms        # check for first call... memory allocation
        goflag = 1
        if(self.errorflag == 1):  return
            
        
      # -------------------- read stuff
        while goflag:

           # read a line
            line = self.fid.readline()
         
           # ------------------------ check for end of file??? 
            if not line:
                self.file_close()

              # open a new file, if possible
                if(self.filenum < len(obj.filename)-1):
                   
                    self.file_open()

                    if(self.errorflag == 1): 
                        return
                    continue
                else:
                    print( "-------------------------------------------")
                    print( "NOTE/WARNING/ERROR:                        ")
                    print( "        Timestep is not contained in user   ")
                    print( "        specified file list of files. This  ")
                    print( "        could be end of files and hence data.")
                    print( "        Exiting softly         ------------")
                    print( "-------------------------------------------")                    
                    self.errorflag = 1;
                    break


           # ------------------------ start parsing data
            if line.find("ITEM: TIMESTEP") > -1:
                line = self.fid.readline()
                temp = str.split(line)
                self.timestep = int(temp[0])
                print("TIMESTEP:  %d" % (self.timestep) )
                if(correct_timestep > -1):
                    if(self.timestep == tstep):
                         correct_timestep = 1
                    elif(self.timestep > tstep):
                         print("-----------------------------------------")
                         print("ERROR:  input dump file does not contain ")
                         print("         the user defined timestep.      ")
                         print("-----------------------------------------")
                         #obj = obj.delete()
                         obj.errorflag = 1
                         return
                    else:
                         correct_timestep = 0;
                #print (" correct timestep %d | %d  %d " % (correct_timestep, obj.timestep, tstep))
                continue
               

            if line.find("ITEM: NUMBER OF ATOMS") > -1:
                line = self.fid.readline()
                temp = str.split(line)
                self.nratoms = int(temp[0])
                continue

            if line.find("ITEM: BOX BOUNDS") > -1:
                temp = str.split(line)
                self.periodicity[0] = temp[-1] 
                self.periodicity[1] = temp[-2]
                self.periodicity[2] = temp[-3]
                
                if line.find(' xy ') > -1:   #this is triclinic

                    line = self.fid.readline()
                    temp = str.split(line)
                    xl = float(temp[0])
                    xh = float(temp[1])
                    xy = float(temp[2])
                    line = self.fid.readline()
                    temp = str.split(line)
                    yl = float(temp[0])
                    yh = float(temp[1])
                    xz = float(temp[2])
                    line = self.fid.readline()
                    temp = str.split(line)
                    zl = float(temp[0])
                    zh = float(temp[1])
                    yz = float(temp[2])
                    self.box = [xl,xh,yl,yh,zl,zh,xy,xz,yz]

                else:

                    line = self.fid.readline()
                    temp = str.split(line)
                    xl = float(temp[0])
                    xh = float(temp[1])
                    line = self.fid.readline()
                    temp = str.split(line)
                    yl = float(temp[0])
                    yh = float(temp[1])
                    line = self.fid.readline()
                    temp = str.split(line)
                    zl = float(temp[0])
                    zh = float(temp[1])
                    self.box = [xl,xh,yl,yh,zl,zh]


                continue


            if line.find("ITEM: ATOMS") > -1:
               
                 temp = str.split(line)
                 ll = temp[2:]
                 N = len(ll)
                 who = np.zeros(shape=(N,1))
                 aid = self.ismember(ll, 'id')

                 if aid == -1:
                     print("-----------------------------------------")
                     print("WARNING: input dump file does not contain ")
                     print("         an atom id number. Atoms cannot ")
                     print("         be sorted. This is an issue when ")
                     print("         comparing two timesteps. ")
                     print("-----------------------------------------")
                 else:
                     who[aid] = 1
                            
                 if (correct_timestep <= 0 and tstep != -1):  # completely skip everything else
                     for ln in range(self.nratoms):
                         tmp = self.fid.readline()
                     #print "skipping"
                     continue

                 self.fields.append('id')       
                #  these are the items it can currently support
                 xid = self.ismember(ll,'x')
                 if xid != -1: self.fields.append('x')
                 yid = self.ismember(ll,'y')
                 if yid != -1: self.fields.append('y')
                 zid = self.ismember(ll,'z')
                 if zid != -1: self.fields.append('z')
                 xsid = self.ismember(ll,'xs')
                 if xsid != -1: self.fields.append('xs')
                 ysid = self.ismember(ll,'ys')
                 if ysid != -1: self.fields.append('ys')
                 zsid = self.ismember(ll,'zs')
                 if zsid != -1: self.fields.append('zs')
                 vxid = self.ismember(ll,'vx')
                 if vxid != -1: self.fields.append('vx')
                 vyid = self.ismember(ll,'vy')
                 if vyid != -1: self.fields.append('vy')
                 vzid = self.ismember(ll,'vz')
                 if vzid != -1: self.fields.append('vz')
                 tid = self.ismember(ll,'type')
                 if tid != -1: self.fields.append('type')
                 qid = self.ismember(ll,'q')
                 #print('tester--------- ',qid)
                 if qid != -1: self.fields.append('q')
                 mid = self.ismember(ll,'mol')
                 if mid != -1: self.fields.append('mol')


                 if xid > -1: 
                     who[xid] = 1
                     present.append('xid')
                 if yid > -1: 
                     who[yid] = 1
                     present.append('yid')
                 if zid > -1: 
                     who[zid] = 1
                     present.append('zid')
                 if xsid > -1: 
                     who[xsid] = 1
                     present.append('xsid')
                 if ysid > -1: 
                     who[ysid] = 1
                     present.append('ysid')
                 if zsid > -1: 
                     present.append('zsid')
                     who[zsid] = 1
                 if vxid > -1: 
                     who[vxid] = 1
                     present.append('vxid')
                 if vyid > -1: 
                     who[vyid] = 1
                     present.append('vyid')
                 if vzid > -1: 
                     who[vzid] = 1
                     present.append('vzid')
                 if tid > -1: 
                     who[tid] = 1
                     present.append('tid')
                 if qid > -1: 
                     who[qid] = 1
                     present.append('qid')
                 if mid > -1: 
                     who[mid] = 1
                     present.append('mid')

#                 print (aid, xid, yid, zid, xsid, ysid, zsid, vxid, vyid, vzid, tid, qid, mid)
#                 print present


              # ---------------- add non-standard per atoms fields
                 unused_ids = []
                 for jj, val in enumerate(who):
                      if val == 0:
                           unused_ids.append(jj)

                 
                 if(ttest == -1):
                      if(len(unused_ids) > 0):
                           self.extras = np.empty(shape=(self.nratoms,len(unused_ids)))
                           for jj, val in enumerate(unused_ids):
                                 self.extra_fields.append(ll[unused_ids[jj]])

             # do we need to reallocate arrays?
                 if self.nratoms != nratoms_old:
                     self.mem_allocate(present)     # custom will need to be included here ...
                     

              # read atoms
                 for ln in range(self.nratoms):   
                     line = self.fid.readline()
                     tmp = str.split(line)
                     
                     if aid > -1:
                         atomid = int(tmp[aid])-1                         
                     else:
                         atomid = ln 

                     if xid > -1:
                         self.coordinate[atomid,0] = float(tmp[xid])
                         self.coordinate[atomid,1] = float(tmp[yid])
                         self.coordinate[atomid,2] = float(tmp[zid])
                      
                     if xsid > -1:
                         #obj.coordinate[atomid,0] = float(tmp[xsid])*(obj.box[1]-obj.box[0])+obj.box[0]
                         #obj.coordinate[atomid,1] = float(tmp[ysid])*(obj.box[3]-obj.box[2])+obj.box[2]
                         #obj.coordinate[atomid,2] = float(tmp[zsid])*(obj.box[5]-obj.box[4])+obj.box[4]
                        # ------------- leave them scaled
                         self.coordinate[atomid,0] = float(tmp[xsid])
                         self.coordinate[atomid,1] = float(tmp[ysid])
                         self.coordinate[atomid,2] = float(tmp[zsid])
                           
                     if vxid > -1:
                         self.velocity[atomid,0] = float(tmp[vxid])
                         self.velocity[atomid,1] = float(tmp[vyid])
                         self.velocity[atomid,2] = float(tmp[vzid])
                            
                     if tid > -1:
                         self.type[atomid] = int(tmp[tid])
                            
                     if qid > -1:
                         self.charge[atomid] = float(tmp[qid])
                         #print(obj.charge[atomid], tmp[qid], float(tmp[qid]))                            


                     if mid > -1:
                         self.molecule[atomid] = int(tmp[mid])
                            
                     ex = max([aid, xid+2, xsid+2, vxid+2, tid, qid, mid])+1

                     if len(unused_ids) > 0:
                         for mm, val in enumerate(unused_ids):
                             self.extras[atomid,mm] = float(tmp[mm+ex])
          
            self.nrtypes = int(max(self.type))
            goflag = 0

        
        
        
        
        return self
    
    
    
    
    def load_data(self,tstep):
        # initialize the first step("dump") ... or load data ("lammps")
        # lammps input data file
        # All arrays are returned as sorted by "id" number        

        if self.inputFileType == "lammps":
            lines = self.fid.readlines()
            ## PARSE EACH LINE FOR REQUIRED INFO
            for i in range(0, len(lines)):
                line = lines[i]
                if 'atoms' in line:
                    nratoms = [x.strip() for x in line.split()] #split the line and remove white space
                    self.nratoms = int(nratoms[0])
                if 'types' in line:
                    nrts = [x.strip() for x in line.split()]
                    self.nrtypes = int(nrts[0])
                    self.masses = np.zeros((self.nrtypes, 1))
                if 'bonds' in line:
                    nratoms = [x.strip() for x in line.split()] 
                    self.nrbonds = float(nratoms[0])
                if 'angles' in line:
                    nratoms = [x.strip() for x in line.split()] 
                    self.nrangles = float(nratoms[0])
                if 'dihedrals' in line:
                    nratoms = [x.strip() for x in line.split()]
                    self.nrdihedrals = float(nratoms[0])
                if 'impropers' in line:
                    nratoms = [x.strip() for x in line.split()]
                    self.nrimpropers = float(nratoms[0])
                if 'xlo xhi' in line:
                    ln = [x.strip() for x in line.split()]
                    xl = float(ln[0])
                    xh = float(ln[1])
                if 'ylo yhi' in line:
                    ln = [x.strip() for x in line.split()]
                    yl = float(ln[0])
                    yh = float(ln[1])
                if 'zlo zhi' in line:
                    ln = [x.strip() for x in line.split()]
                    zl = float(ln[0])
                    zh = float(ln[1])
                if 'Masses' in line:                   
                    i +=1 # ignore the empty line                  
                    jj = 0
                    while jj < self.nrtypes:
                        line = lines[i+1]
                        i +=1 
                        if line:
                            ll = [x.strip() for x in line.split()]
                            self.masses[int(ll[0])-1] = float(ll[1])
                        else:
                            break
                        jj +=1
                ## FURTHER PARSE STORED INFO
                if 'full' in self.atom_style:
                    c_indx = np.arange(5,8) -1
                    t_indx = 2
                    m_indx = 1
                    q_indx = 3
                    f_indx = np.arange(8,11)-1
                    # frmt = '%f %f %f %f %f %f %f %f %f %f'
                if 'molecular' in self.atom_style:
                    c_indx = np.arange(4,7)-1
                    t_indx = 2
                    m_indx = 1
                    q_indx = -1
                    f_indx = np.arange(7,10)-1
                    # frmt = '%f %f %f %f %f %f %f %f %f'
                  
                if 'Atoms' in line:
                    i +=1
                    
                    #initialize arrays, but only the first time
                    # if self.nratoms > 0:
                    self.coordinate = np.zeros((self.nratoms,3));
                    self.type = np.zeros((self.nratoms,1));
                    self.molecule = np.zeros((self.nratoms,1));
                    self.charge = np.zeros((self.nratoms,1));
                    self.image_flag = np.zeros((self.nratoms,3));
                    
                    # determine number of cols
                    line = lines[i+1]
                    i +=1 
                    ll = [x.strip() for x in line.split()]
                    ll = [float(x) for x in ll]
                  
                    N = len(np.where(~ np.isnan(ll))[0])
                    # frmt = np.tile('%f',(1,N))
                    # f2 = np.tile('%s',(1,len(ll)-N))
                    # frmt = np.concatenate((frmt,f2), axis = 1)
                    
                    
                    # matlab obj.coordinate(str2double(ll(1)),:) = str2double(ll(c_indx));
                    # python indexing starts from 0
                    self.coordinate[int(ll[0])-1,:] = ll[c_indx[0]:c_indx[-1]+1]
                    self.type[int(ll[0]-1)]= ll[t_indx]
                    self.molecule[int(ll[0]-1)] = ll[m_indx]
                    if q_indx != -1:
                        self.charge[int(ll[0]-1)]= ll[q_indx]
                    if N > 7:
                        self.image_flag[int(ll[0])-1,:] = ll[f_indx[0]:f_indx[-1]+1]
                    
                    txt = np.genfromtxt(lines[i+1:]) # get the remaining lines
                    ids = np.array([int(txt[j][0]) for j in range(0,len(txt))])
                    
                    # c1 = np.array([txt[j][c_indx[0]] for j in range(0,len(txt))])
                    # c2 = np.array([txt[j][c_indx[1]] for j in range(0,len(txt))])
                    # c3 = np.array([txt[j][c_indx[2]] for j in range(0,len(txt))])
                    # self.coordinate[ids - 1,:] = np.array[c1, c2, c3]
                    
                    for j in range(0,len(ids)):
                        self.coordinate[ids[j]-1] = [ txt[j][c_indx[0]], txt[j][c_indx[1]], txt[j][c_indx[2]]]
        
                if 'Velocities' in line:
                    self.velocity = np.zeros(shape=[self.nratoms,3])
                    self.fid.readline()
                    
                    for ll in range(0,self.nratoms):
                        temp = str.split(self.fid.readline())
                        idx = int(temp[0])-1
                        self.velocity[idx,:] = [float(temp[1]), float(temp[2]), float(temp[3])]
        
            # except IOError:
            #     sys.stderr.write( "lmp_class.py:read_lammps - Error: Could not open %s\n", (self.filename) )
            #     sys.exit(-1)
            # return
        
        return self
    
    
    def ismember(self,lst, string):
        rtn = -1         
        ii = 0
        while ii < len(lst):
            if string == lst[ii]:
                rtn = ii
                return rtn
            ii += 1
        return rtn
    
    def mem_allocate(self,lst):

        #obj.delete()

        # print ("allocating memory")
        N = self.nratoms
        self.coordinate = np.zeros(shape=(N,3))
        if self.ismember(lst,'vxid') > -1: self.velocity = np.zeros(shape=(N,3))
        if self.ismember(lst,'mid') > -1: self.molecule = np.zeros(shape=(N,1), dtype=int)
        if self.ismember(lst,'qid') > -1: self.charge = np.zeros(shape=(N,1), dtype=float)
        if self.ismember(lst,'tid') > -1: self.type = np.zeros(shape=(N,1), dtype=np.int8)
        if self.ismember(lst,'ifid') > -1: self.image_flag = np.zeros(shape=(N,3), dtype=np.int8)

        return
    
    def open_file(self):
        
        print ("Opening ", self.filename, "\n")
        
        try:
            fid = open(self.filename,'r')
        
        except OSError:
            print ("Could not open/read file:", self.filename)
            sys.exit()
          
        return fid
    
    def close_file(self):
        self.fid.close()