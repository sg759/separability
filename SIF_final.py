
# Python implementation completed by Swati Gupta in September 2022.

# Geometrical method, creating toy config, MD application code adapted from
# the Matlab implementation by Dr. Mark Wilson and Dr. Scott J. Grutzik (2019)

# This code builds on some previous work done in Matlab by Dr. Derek Warner and Grant West.
#
# Uses PSIFA (Pattern Search and Implicit Filtering Algorithm), a derivative-free 
# optimization algorithm developed by Ferreira, Ehrhardt and Santos (2017).

from inputModulesHelpers import *

'''
This module calculates the cracktip and the associated K-field values using
three approaches:
    - Separability approach
    - Geometrical or Last Surviving Bond (LSB) 
    - Displacement Correlation (DC)
     
INPUTS:
   application: synData = 1 -> Toy Model that generates synthetic data
                synData = 0  -> reference coordinates and displacement field provided
                
   r2:          outer radius of annulus
   
   alpha:       scaling factor that determines the size of the inner radius of the annulus
                
   guess:       initial guess of crack top location. Should be an array of two
                coordinates, eg, [x,y]
   
   constants:   values of elastic constants. The default values are (in order)
                poisson =  0.25          # poisson's ratio
                E = 70                   # Young's modulus
                shear_modulus = 0.25  
                
                These are used to calculate the associated K-field values.
  
   For Toy Model (Synthetic Data):   
          
   domainSize:  #points for synthetic data in Toy Model. Can be a scalar or
                an array of #points
                
    crack:      cracktip location used in synthetic data. Should be an array of two
                coordinates, eg, [x,y]
                
    shape:      0: on lattice - single radius
                1: random in a box
                2:ring
                3: on lattice
    
    noise:      magnitude of noise in the synthetic data in Toy Model
    K_app:      applied K field
     
   
   For MD/DIC/FE Data:
   coords_ref:  reference coordinates
   coords:      coordinates after displacement

Control variables: 
   geo: use Geometrical method (geo = 1) or not (geo = 0)
   Sep: use Separability method (Sep = 1) or not (Sep = 0)
   DC: use Displacement Corr method (DC = 1) or not (DC = 0)
   PS: use Pattern Search (PS = 1) or not (PS = 0)

OUTPUTS:
   cracktip1:    crack tip location found using geometrical method
   K_field1:     K field associated with cracktip found using geometrical method

   cracktip2:   crack tip location found using seperability 
   K_field2:    K field associated with the cracktip found using seperability 
   
   cracktip3:   crack tip location found using DC 
   K_field3:    K field associated with the cracktip found using DC 

'''
def SIF_projection(synData = 1, shape = 3, domainSize=10000,r2 = 50,alpha=1.5,noise = 0, 
                   K_app = [0.8, 0, 0], coords = [], coords_ref =[],config = 0,guess = [10, 10], 
                   crack = [0,0], h=5,geo = 0,Sep = 1, DC = 0, PS = 1, constants = []):
            
    global R1, R2, boxsize, tip_guess, poisson, E, kolosov, shear_modulus
    
    R2 = r2
    R1 = r2/(alpha**2)
    boxsize = 4*R2 
    # pdb.set_trace()
    tip_guess = guess
    K_app = np.array(K_app)
    
    plot_configs = 0
    
    if max(np.shape(constants)) != 0:
        poisson = constants[0]
        E = constants[1] 
        shear_modulus = constants[2]
        kolosov = 3-4*poisson
    
    if (synData == 1):
        
        print("\n\n Toy Model \n")
        
        [coords,r,theta,dx,coords_ref,K_field, ids]= ToyModel(synData,shape,alpha,domainSize,noise,crack,h,K_app, geo,Sep,PS,DC)
        lb = np.array([-boxsize/2, -boxsize/2]) #upper bound for pattern search
        ub = np.array([boxsize/2, boxsize/2]) #lower bound for pattern search        
    else: 
        print("\n\n Displacement data provided \n")
        # pdb.set_trace()
        [theta,r] = cart2pol(coords[:,0]-guess[0], coords[:,1]-guess[1]);               
        theta = theta + pi/2;
        theta[theta > pi] =  theta[theta > pi] - 2*pi;
        
        dx = dx = np.array(coords - coords_ref)
        if synData == 2:
            lb = np.array([min(coords_ref[:,0]), 20+(min(coords_ref[:,1])/4)])
            ub = np.array([max(coords_ref[:,0]/2), 30+(max(coords_ref[:,1])/4)])
        else:
            lb = np.array([min(coords[:,0]), min(coords[:,1])])
            ub = np.array([max(coords[:,0]), max(coords[:,1])])
        print ("\n guess: ", guess)
        K_app = []
    
    # initilize variables
    cracktip1 = [-100,100] # Geo
    cracktip2 =[-100,100] # Sep
    cracktip3 = [-100,100] # DC
    K_field1 = np.zeros((1,6))  # Geo
    K_field2 = np.zeros((1,6))  # Sep
    K_field3 = np.zeros((1,6))  # DC
    
    if geo:
        # find the cracktip using last surviving bond (geometric method)
        print('CALCULATING GEOMETRIC CRACK TIP LOCATION\n')
        edges = find_boundary(np.array([coords[:,0], coords[:,1]])); 
       
        B0 = [x[0] for x in edges]; # x coords
        B1 = [x[1] for x in edges]; # Y coords
        B = np.vstack([B0,B1]);
    
        [boxx, boxy] = gen_box_poly(boxsize,tip_guess[0],tip_guess[1]); # generate a polygon about the previous crack tip location
        if(plot_configs):
            cracktip1, fig1= locate_cracktip(coords, B, boxx, boxy, tip_guess, plot_configs, synData)
        else:
            cracktip1= locate_cracktip(coords, B, boxx, boxy, tip_guess, plot_configs, synData)
        print('CRACKTIP LOCATION (x,y) geometric: \n', cracktip1)      
        
        K_field1, ids1 = calcKField(coords_ref,dx,cracktip1)
        
    
    if Sep & PS:
        print('CALCULATING SEPARABLE CRACK TIP LOCATION using Pattern Search\n')
        
        emp = np.array([])
        funcSep = lambda x:costFuncSep(x,coords_ref,dx,alpha)
        cracktip2, cfSep, iterations, alpha, evaluations, history, exit_cause = \
            psifa(funcSep,np.array(guess), emp,emp,emp,2,lb,ub)
        
        print('CRACKTIP LOCATION (x,y) separability approach with  Pattern Search: \n', cracktip2)
        K_field2, ids2 = calcKField(coords_ref,dx,cracktip2)
    elif (Sep == 1 and PS == 0):
        
        print('CALCULATING SEPARABLE CRACK TIP LOCATION\n')
        cracktip2 = tip_guess      
        cracktip2,cfSep = find_sep_cracktip(coords_ref,dx,h,alpha,noise,K_app, synData)
       
        cracktip2 = [cracktip2[1], cracktip2[0]]
        print('CRACKTIP LOCATION USING SEPARABLE APPROACH\: \n', cracktip2)
        K_field2, ids2 = calcKField(coords_ref,dx,cracktip2)
        
    if DC:
        K_guess = [1,0.5,0.03]
        cracktip3, min_ori, K_field3,cfDC,evaluations = dispCorr(dx,coords,K_guess,synData)
        K_field3, ids3 = calcKField(coords_ref,dx,cracktip3)
        print('CRACKTIP LOCATION USING Displacement Correlation(DC): \n', cracktip3)
        
    if(plot_configs and geo):
        plt.figure(fig1)
        plt.plot(cracktip2[0], cracktip2[1],marker = "P", markersize = 7)
        label = f"Sep({cracktip2[0]:0.3f}, {cracktip2[1]:0.3f})" 
        plt.annotate(label,(cracktip2[0], cracktip2[1]),textcoords="offset points",xytext=(0,20),ha='left')
        # plt.grid(True)
        plt.title("Crack Tip Locations for Toy Model, N =" + str( domainSize))
        plt.show()
    elif (plot_configs and DC):
        plt.figure()
        plt.plot(coords[:,0],coords[:,1])
        plt.plot(cracktip3[0], cracktip3[1],marker = "P", markersize = 7)
        label = f"DC({cracktip3[0]:0.3f}, {cracktip3[1]:0.3f})" 
        plt.annotate(label,(cracktip3[0], cracktip3[1]),textcoords="offset points",xytext=(0,20),ha='left')
        # plt.grid(True)
        plt.title("Crack Tip Location")
        plt.show()

    if synData == 2: #MD
        if geo:
            K_field1 = [K_field1[0], K_field1[1]*10**-2,K_field1[2]*10**-2,K_field1[3]*10**-7,K_field1[4],K_field1[5]]
        if Sep:
            K_field2 = [K_field2[0], K_field2[1]*10**-2,K_field2[2]*10**-2,K_field2[3]*10**-7,K_field2[4],K_field2[5]]
        if DC:
            K_field3 = [K_field3[0], K_field3[1]*10**-2,K_field3[2]*10**-2,K_field3[3]*10**-7,K_field3[4],K_field3[5]]
            

    print('\n---------------------------------------------------------------------------------\n')
    if Sep:
        print('NUMBER OF POINTS SAMPLED  \n', np.size(ids2) );
        print(' K1          K2          T      \n')
        print('K Field Separability:', K_field2[2:5])
    if DC:
        print('K Field DC:', K_field3[2:5])
    if geo:
        print('K Field Geometric:', K_field1[2:5])
    # print ["{0:0.3f}".format(i) for i in K_field] 
    # print('Applied:',K_app)
    # print[ "{0:0.3f}".format(i) for i in K_app]
              
    return cracktip1,cracktip2,cracktip3,K_field1,K_field2,K_field3


def ToyModel(choice,shape,alpha,numPoints,noise,crack,h, K_app, geo,Sep,PS, DC):
    ###############################################################################################
    ### Build a random initial config by apply the known displacement to a random set of points ###
    ###############################################################################################
    
    global  N, boxsize, R_0, nmax, poisson, E, tol, kolosov, shear_modulus,nrpoints
 
    
    N = numPoints
    random_magnitude = noise      # magnitude of random noise  
    cracktip = crack        # where are you locating the crack tip
    plot_configs = 1
    
    print("\n\nBUILDING INITIAL CONFIG\n")#, N[i]);
    
    # assure the displacements are in reference to the correct coordinate system

    coords,r,theta,dx,coords_ref = gen_random_config(shape, K_app, int(N), random_magnitude, cracktip, plot_configs);
   
    print('CALCULATING STRESS FIELD ASSUMING USER SPECIFIED CRACK TIP LOCATION\n')
    K_field, ids = calcKField(coords,dx,cracktip)
    
    return coords,r,theta,dx,coords_ref, K_field, ids 

def calcKField(coords,dx,tip):
    ###########################################################
    ### Compute K Field asscociated with given tip location ###
    ###########################################################
    
    [theta,r] = cart2pol(coords[:,0]-tip[0], coords[:,1]-tip[1]);             
    theta = theta + pi/2;
    theta[theta>pi] =  theta[theta>pi] - 2*pi;      
    # pdb.set_trace()
    [Area, ids] = determine_area3(r,theta,1);  
    
    # build the vector and matrix for the linear system:  Ax = b
    A = build_matrix(ids,r,theta,Area);
    b = build_lhs(ids, dx[:,0], dx[:,1], r, theta, Area)
    # solve it and convert to stress field   
    
    K_field = np.linalg.inv(A.T @ A) @ A.T @ b;
    K_field = np.multiply(K_field.T , np.array([1 , 1,  math.sqrt(2*math.pi),\
                                -math.sqrt(2*math.pi) , 4 , np.ones(np.size(K_field)-5)[0]]))
    K_field = K_field[0]
    
    return K_field, ids


def find_sep_cracktip(coords,DX,h,alpha,noise=0,K_app='unknown',synData = 0):
    ########################################################################
    # SEPARABILITY APPROACH                                                #
    # This function finds the cracktip by finding the minimum value of the #
    # separability cost function, as decribed in our paper                 #
    ########################################################################
    
    
    global R1, R2, random_or_all, nrpoints, boxsize
    
    if synData != 1:
        bbox = [min(coords[:,0]), max(coords[:,0]),20+(min(coords[:,1])/4),(max(coords[:,1])/2)]
         
    else:
        bbox = [-boxsize/2, boxsize/2, -boxsize/2, boxsize/2];
   
    xx = np.arange(bbox[0],bbox[1]+1,h) #crack tip positions to be examined in x direction
    yy = np.arange(bbox[2],bbox[3]+1,h) #crack tip positions to be examined in y direction
    
    num_subannuli = 2; #number of sub-annuli to be used
    
    # extract displacement data
    dx = DX[:,0] 
    dy = DX[:,1]
    
    [XX,YY] = np.meshgrid(yy,xx); # for plotting
    shapeXX = np.shape(XX)
    # Flatten: in order to index a 2D array linearly a la Matlab
    XX = XX.flatten('F')
    YY = YY.flatten('F')
    
    Rinner = [R1, alpha*R1]; # inner radius of sub-annuli
    cf = np.zeros((XX.size,1))

    #loop over all crack tip positions
    for jj in range(0,XX.size):
        
        [theta,r] = cart2pol(coords[:,0]-XX[jj], coords[:,1]-YY[jj]);
        avg_disp = np.zeros((num_subannuli, 2))
        
        #loop over sub annuli
        for kk in range(0,num_subannuli):
            ids = np.where((r >= Rinner[kk]) & (r < alpha*Rinner[kk]))[0]  
            # Note: brackets within the above condition between '&' are important
            if len(ids) == 0:  
                # Make the value nan so that
                avg_disp[kk,0] = np.NAN
                avg_disp[kk,1] = np.NAN 
            else: 
               
                avg_disp[kk,0] = np.mean(dx[ids])#, axis = 0)
                avg_disp[kk,1] = np.mean(dy[ids])#, axis = 0)       
        
        # compute C1 and C2
        C1b = (1/num_subannuli)*np.sum(avg_disp[:,0]/np.sqrt(Rinner))
        C2b = (1/num_subannuli)*np.sum(avg_disp[:,1]/np.sqrt(Rinner)) 
        
        # compute the cost function for current crack tip position
        cf[jj] =  (1/num_subannuli) * np.sum((avg_disp[:,0]-C1b*np.sqrt(Rinner))**2 + (avg_disp[:,1]-C2b*np.sqrt(Rinner))**2,0)
                                                
    mv = min(cf[:])
    val = np.reshape(cf,shapeXX)
    cc,rr = np.where(val == mv)  
    if synData == 1:
        cracktip = [-yy[rr[0]], xx[cc[0]]]   
    else:
        cracktip =[xx[cc[0]], yy[rr[0]] ]  
    
    return cracktip,cf

def costFuncSep(guess,coords,DX,alpha):
    ########################################################################
    # Objective function for minimizng SEPARABILITY cost function 
    ########################################################################
    
    global R1, random_or_all, nrpoints
    
    num_subannuli = 2; #number of sub-annuli to be used
    
    # extract displacement data
    dx = DX[:,0] 
    dy = DX[:,1]
    
    avg_disp = np.zeros((num_subannuli, 2)) #average displacement components for each subannuli
    Rinner = [R1, alpha*R1]; # inner radius of sub-annuli
       
    [theta,r] = cart2pol(coords[:,0]-guess[0], coords[:,1]-guess[1]);
    
    for kk in range(0,num_subannuli):
        ids = np.where((r >= Rinner[kk]) & (r < alpha*Rinner[kk]))[0]  
        # Note: brackets within the above condition between '&' are important
       
        avg_disp[kk,0] = np.mean(dx[ids])
        avg_disp[kk,1] = np.mean(dy[ids])     
    
    # compute C1 and C2
    C1b = (1/num_subannuli)*np.sum(avg_disp[:,0]/np.sqrt(Rinner))
    C2b = (1/num_subannuli)*np.sum(avg_disp[:,1]/np.sqrt(Rinner)) 
    
    # compute the cost function for current crack tip position
    cf =  (1/num_subannuli) * np.sum((avg_disp[:,0]-C1b*np.sqrt(Rinner))**2 + (avg_disp[:,1]-C2b*np.sqrt(Rinner))**2,0)
    eval_success = True
    return cf,eval_success

def dispCorr(dx,coords,K_guess,app):
    
    emp = np.array([])
    
    ori = 0;
    K1 = K_guess[0];
    K2 = K_guess[1];
    T = K_guess[2];
    
    params0 = np.array([ori, tip_guess[0], tip_guess[1], K1,K2,T])
    
    ## DIC bounds
    if app ==0:
        lb = np.array([-pi/2, -boxsize/2, -boxsize/2, 0, -5, -5]) #lower bounds of params
        ub = np.array([pi/2, boxsize/2, boxsize/2, 2, 2, 1]) # upper bounds
    else:
        lb = np.array([-pi/2, min(coords[:,0]), 10+min(coords[:,1]), 0, -5, -5])
        ub = np.array([pi/2,max(coords[:,0]), max(coords[:,1]),2, 2, 1])
    n = 6 #number of minimization problems
    
    func = lambda params:costFuncDC(params,coords,dx)
    x, cfDC, iterations, alpha, evaluations, history, exit_cause = \
                psifa(func,params0, emp,emp,emp,n,lb, ub)
    cracktip = x[1:3]
    min_ori = x[0]
    min_KField = x[3:6]
    
    
    return cracktip, min_ori, min_KField,cfDC,evaluations

def costFuncDC(params,coords, dx):
    ########################################################################
    # Objective function for minimizing orientation, crack tip guess and
    # Kfield = [K1,K2,T] for Displacement Correlation  cost function
    ########################################################################
    # parse input arguments
    orientation = params[0]
    guess = params[1:3]
    K_field = params[3:6]
    # pdb.set_trace()
    [theta,r] = cart2pol(coords[:,0]-guess[0], coords[:,1]-guess[1])
    
    # find points within the annulus
    ids = np.where((r >= R1) & (r < R2))[0]  
    if (max(ids.shape) == 0):
        cf = np.NaN;
        return cf,False
  
    if(random_or_all == 'r'):
        rn = randperm(max(ids.shape),nrpoints)
        ids = ids[rn]

    # apply given crack orientation
    theta_rot_test = theta + pi/2 + orientation; 
    theta_rot_test[theta_rot_test > pi] =  theta_rot_test[theta_rot_test > pi] - 2*pi
    # get Williams Expansion functions
    dx1 = u_exp(theta_rot_test,r,1,1).T
    dx2 = u_exp(theta_rot_test,r,1,2).T
    dx3 = u_exp(theta_rot_test,r,2,1).T
    
    
    # multiply with guess for Kfield to get William's Coefficients
    A_app = K_field / np.array([np.sqrt(2*pi), -np.sqrt(2*pi), 4])

    dx_rot = A_app[0] * np.array([dx1[:,1], -dx1[:,0]]) + A_app[1] * np.array([dx2[:,1], -dx2[:,0]]) \
             + A_app[2] * np.array([dx3[:,1], -dx3[:,0]])
             
    dx_rot = dx_rot.T

    #calculate residual
    cf = (1/max(ids.shape)) * sum((dx_rot[ids,1] - dx[ids,1])**2 + (dx_rot[ids,0] - dx[ids,0])**2)
    eval_success = True
    return cf,eval_success
    

def determine_area3(r,t,fignum):
    ########################################################################
    # Function returns average area per point (scalar) and the ids of those
    # contained in the annulus
    ########################################################################
    
    ids = np.where((r <= R2) & (r > R1 )& (t <= tmax) & (t > tmin))[0];
    
    lengthIDs = max(ids.shape);
    
    # fix the number of points selected 
    if(random_or_all == 'r'):
        if(lengthIDs < nrpoints): #length(ids) = max(ids.shape)
            print('%length(ids)  nrpoints \n', lengthIDs, nrpoints);
        rn = randperm(lengthIDs,nrpoints);
        #rn = np.random.permutation(length,nrpoints);
        ids = ids[rn];
        
    newLength = max(ids.shape)
    #plot([0 R2.*cos(pi/2 + 0.05)], [0 R2.*sin(pi/2 + 0.05)], 'linestyle','-','color','k','linewidth',2)
    #A = 0.5*(R2^2 - R1^2)*(tmax-0.05-(tmin+0.05));
    A = 0.5*(R2**2 - R1**2)*(tmax-tmin);
  
    AA = A/newLength;
  
    return AA, ids;
     

def locate_cracktip(coords, B, boxx, boxy, cracktip, plot_stuff,synData = 1):
    ############################################################################
    # Last Surviving Bond or GEOMETRIC approach
    # function  finds the lowest magnitude x coordinate  
    ############################################################################
  
    n = len(boxx)
    # coordinates of B outline the domain including the shape
    # of the crack, in otherwords B is index of all surface atoms
    B=B[0]
    
    # Find all surface atoms within a box near the crack tip
    IN, ON = inpolygon(coords[B,0], coords[B,1], boxx[0:n-1],boxy[0:n-1]);
   
    idIN = np.where(IN);
    idON = np.where(ON);
    
    bp = idIN[0];
     
    if( bp.any() ):
        # of those inside the box, find the smallest in the y direction
        zz = np.array([coords[B[bp],0], coords[B[bp],1]])
        l = max(np.shape(zz))
        z = [(zz[0][i],zz[1][i]) for i in range(0,l)]
        if synData == 1:
            z.sort(key = lambda z: z[1])
        else:
            z.sort(key = lambda z: z[0])
        # average value of coords with 3 smallest y values
        xy = np.mean(z[0:3],axis=0) 
        if synData == 1:
            cracktip = [-xy[1], xy[0]]
        else:
            cracktip = [xy[0], xy[1]] 
    
    
    if(plot_stuff):
        fig1 = plt.figure()
        if synData == 1:
            plt.plot(-coords[B,1], coords[B,0],'go')
            plt.plot(-coords[B[bp],1], coords[B[bp],0],'ro')
        else:
            plt.plot(coords[B,0], coords[B,1],'go')
            plt.plot(coords[B[bp],0], coords[B[bp],1],'ro')
        plt.plot(cracktip[0], cracktip[1],marker = "X", markersize = 9)
        label = f"LSB({cracktip[0]:0.3f}, {cracktip[1]:0.3f})" 
        plt.annotate(label,(cracktip[0], cracktip[1]),textcoords="offset points",xytext=(0,-15),ha='left')
        plt.plot(boxx, boxy)#'color',[0.4, 0.4, 0.4],'linestyle','--', 'linewidth',2)
        plt.xlabel('X direction') 
        plt.ylabel('Y direction') 
        # plt.title("Last Surving Bond Crack Tip Location (PYTHON) for N = 30,000")
        plt.show()
        return cracktip, fig1
    else:
        return cracktip

def build_matrix(ids,r,theta,A):
############################################################################
# Function will build the matrix of the linear system, defined as
# "S", then overlap matrix in equation 11 of Wilson et al. "Continuum stress intensity
# factors from atomistic fracture simulations", CMAME(2019)
# 
# This function is used in both toy-model and MD analysis
############################################################################   
    mat = np.zeros((nmax*2+2,nmax*2+2));
    
    N = np.shape(mat)[1];
    # In python, range returns a sequence of numbers, in increments of 1 (by default), 
    # and stops before N+1, thus including N.
    for rr in range(1,N+1):
        for cc in range(1,N+1):
            n1 = math.ceil(cc/2)-1;            
            n2 = math.ceil(rr/2)-1;           
            kk1 = 2 - (cc % 2) #mod(cc,2);            
            kk2 = 2 -(rr % 2);            
            dx1_ = u_exp(theta[ids],r[ids],n1,kk1);          
            dx2_ = u_exp(theta[ids],r[ids],n2,kk2);
            
            # Convert dx1 and dx2 to 2D array
            dx1 = (np.reshape(dx1_, (2,-1))).T;            
            dx2 = (np.reshape(dx2_, (2,-1))).T;
                         
            ## For array, ``*`` means element-wise multiplication, while ``@`` means matrix multiplication
            ## np.einsum('ij,ij->i',A, B) treats the rows of A and B as vectors and 
            ## returns the dot products of corresponding rows. i.e same as
            ## Matlab's dot(A,B,2). Note that Python's np.dot() does not have an
            ## option to specify the dimension across which the function should
            ## operate along.
     
            mat[rr-1,cc-1] = mat[rr-1,cc-1] + np.sum(np.einsum('ij,ij->i',dx1, dx2) * A);
   
    return mat;

def build_lhs(ids, dy, dx, r, theta, A, choice=1):
    #####################################################################
    # Function will build the vector of the linear system b = Ax, defined
    # as "b" the overlap matrix in equation 11 of Wilson et al. "Continuum 
    # stress intensity factors from atomistic fracture simulations", CMAME(2019)
    #
    # function used in both toy-model and MD analysis
    #
    # NOTE that x_sim -> y_LEFM and y_sim -> -x_LEFM
    # hence, dx and dy variables are switched 
    #####################################################################
    dx = -dx;
    N_ids = np.size(ids);
    br = np.zeros(2 * N_ids,);
   
    # make every other element = value from dx
    if choice == 0:
        br[0: 2 * N_ids: 2] = dx[0][ids]
        br[1: 2 * N_ids: 2] = dy[0][ids]
    else:
        br[0: 2 * N_ids: 2] = dx[ids]
        # make every other element,starting from the second element = value from dy
        br[1: 2 * N_ids: 2] = dy[ids]
         
    N = max(np.shape(br));
    ket = np.zeros((N,1));
    vect = np.zeros((nmax*2+2,1));
    for jj in range(1,nmax*2+2 + 1):
        ii = 0
        kk1 = 2-(jj % 2);
        n1 = math.ceil(jj/2)-1;
   
        for kk in range(0,N,2) :
            ddx = u_exp(theta[ids[ii]],r[ids[ii]],n1,kk1);
            ket[kk] = ddx[0]*A;
            ket[kk+1] = ddx[1]*A;
            ii = ii + 1;
       
        vect[jj-1] = np.dot(br, ket);   
   
    return vect;

def gen_random_config(init_num, K_app, N, random_magnitude, tip, plot_stuff):
    #####################################################################
    # this function will generate various types of initial conditions (toy models)
    # based on user defined variable init_num
    #
    # returns:
    #    coords - matrix of coordinates (N,2)
    #    r -      vector of radial values for each point relative to cracktip
    #    theta -  vector of angles for each point relative to crack face plane
    #    dx -     matrix (N,2) displacment magnitude [dx, dy] for each point, relative to initial conditions 
    
    # function used in both toy-model only
    ##################################################################### 
    
    if (init_num == 1): # random box
    
        coords = 500.*(np.random.rand(N,2)-0.5)
        [theta,r] = cart2pol(coords[:,0]-tip[1], coords[:,1]+tip[0])               
        theta = theta + pi/2
        theta[theta > pi] =  theta[theta > pi] - 2*pi
        
    elif (init_num == 2):  #ring
        
        theta = (np.linspace(-pi,pi,num = N))      
        r = 100.* np.ones((N))
        coords = np.zeros((N,2))    
        [coords[:,0], coords[:,1]] = pol2cart(theta, r)  
        theta = theta + pi/2
        theta[theta > pi] =  theta[theta > pi] - 2*pi
        
    elif (init_num == 3):  #lattice
        XX = 500
        YY = 500
        n = math.floor(np.sqrt(N))
        N = n*n;
        x = np.linspace(-XX/2,XX/2,n)
        y = np.linspace(-YY/2,YY/2,n)
        [X,Y] = np.meshgrid(x,y)
        coords = np.zeros((N,2))    
        coords[:,0] = np.reshape(X,(N), order = 'F')
        coords[:,1] = np.reshape(Y,[N], order = 'F')
        [theta,r] = cart2pol(coords[:,0]-tip[1], coords[:,1]+tip[0])
        theta = theta + pi/2;
        theta[theta > pi] =  theta[theta > pi] - 2*pi
        
    elif (init_num == 0): #lattice - single radius
        XX = 500
        YY = 500
        n = math.floor(np.sqrt(N))
        N = n*n
        x = np.linspace(-XX/2,XX/2,n)
        y = np.linspace(-YY/2,YY/2,n)
        [X,Y] = np.meshgrid(x,y)
        coords = np.zeros((N,2))
        coords[:,0] = np.reshape(X,(N), order = 'F')
        coords[:,1] = np.reshape(Y,(N), order = 'F')
        [theta,r] = cart2pol(coords[:,0]-tip[1], coords[:,1]+tip[0]) 
        theta = theta + pi/2;
        theta[theta > pi] =  theta[theta > pi] - 2*pi  
        
        r = 150.* np.ones((max(theta.shape),1))
        
    else:
        print("\n Invalid Initial Condition \n")
    
    coords_ref = coords;
    A_app = K_app / [np.sqrt(2*pi), -np.sqrt(2*pi), 4]
    dx1 = u_exp(theta,r,1,1)
    dx2 = u_exp(theta,r,1,2)
    dx3 = u_exp(theta,r,2,1)
    rn = random_magnitude *(np.random.rand(N,2)-0.5)
    # When using Python 3.5+, use * for elementwise multiplication 
    # and @ for matrix multiplication
    dx_x = A_app[0] * np.array([dx1[1,:], -dx1[0,:]]) + A_app[1] * np.array([dx2[1,:], -dx2[0,:]]) + A_app[2] * np.array([dx3[1,:], -dx3[0,:]])
  
    ## dx_x is a 3d array. Convert to 2D array, same as rn. 
    ## If A.shape = (a,b,c), D = A.reshape(a,-1) gives D.shape = (a, b*c) 
    dx_xx = np.reshape(dx_x, (2,-1))
    dx = dx_xx.T + rn         

    coords = coords_ref + dx  
 
    if(plot_stuff):
        dd = np.sqrt(dx[:,0]**2 + dx[:,1]**2)
        # plot initial conditions
        
        f1 = plt.figure(1); plt.clf()
        plt.gca().set_aspect('equal', adjustable='box'); # Matlab code: axis equal; or plt.axis('square')
        
        plt.scatter(-coords[:,1], coords[:,0],10 * np.ones((N,1)),dd,marker='o')
        plt.ylabel('Y direction ',fontsize=20)
        plt.xlabel('X direction',fontsize=20)

        h2 = plt.colorbar()
        h2.ax.set_ylabel('Displacement',fontsize=20)
        plt.show()
          
    return coords, r, theta, dx, coords_ref

def u_exp(theta,r,n,K):
    ######################################################################
    # Calculate a displacement field according to the Williams expansion,
    # equations 1-4 of Wilson et al. "Continuum stress intensity
    # factors from atomistic fracture simulations", CMAME(2019)   
    # The returned field is evaluated at points corresponding to r and theta 
    # % Ayatollahi, M. R.  Optics And Lasers in Engineering 90 (2017) 26-33
    # Malikova L.  Fatigue and Fracture of Engineering Materials and
    #               Structures 2015, 38, 91-103
    # Ayatollahi, M. R.  Fatigue and Fracture of Engineering Materials and
    #               Structures 2011, 34, 159-176
    
    # function used in both toy-model and MD analysis
    
    
    # NOTE:   n ranges from 1 to ...N   not[0,N]
    # u - x displacment
    # v - y displacement
    #####################################################################
      
    # converting if(isvector(r)) to python as follows:
    if isinstance(r, (list, np.ndarray, np.matrix)):
        if(~iscolumn(r)):
            r = r.T; 
        if(~iscolumn(theta)): 
            theta = theta.T; 
    else:
        #np.reshape(X,(N), order = 'F');
        r = np.reshape(r, (np.prod(np.size(r)),1),order = "F");
        theta = np.reshape(theta, (np.prod(np.size(theta)),1), order = 'F');

    
    if(K == 1):       # for determination of K1_n values
          u =  (kolosov+(n/2)+(-1)**n) * np.cos( (n/2) * theta ) - (n/2) * np.cos( ((n/2)-2) *theta );
          v =  (kolosov-(n/2)-(-1)**n) * np.sin( (n/2) * theta ) + (n/2) * np.sin( ((n/2)-2) *theta );
        
    elif(K == 2):    # for determination of K2_n values
          u =  (-kolosov-(n/2)+(-1)**n) * np.sin( (n/2) * theta ) + (n/2) * np.sin( ((n/2)-2) * theta );
          v =  ( kolosov-(n/2)+(-1)**n) * np.cos( (n/2) * theta ) + (n/2) * np.cos( ((n/2)-2) * theta );
    
    
    u = u * 1/(2*shear_modulus) * r**(n/2);
    v = v * 1/(2*shear_modulus) * r**(n/2);
    
    displace = np.array([u,v]);
    
    return displace;

