import numpy as np
import matplotlib.pyplot as plt
import Fortran_Code

# Compute the flux terms
def fv_sw_flux_zhang(sim):

    # Some storage fields
    the_flux = np.zeros(sim.sol.shape)
    hs = np.zeros((sim.Nx,sim.Ny,sim.Nz+1))

    # Difference in layer deformations gives layer depths
    hs[:,:,:sim.Nz] = sim.sol[sim.Ih,:,:,:sim.Nz] - sim.sol[sim.Ih,:,:,1:]

    # Loop through each layer and compute flux
    for ii in range(sim.Nz):
      
        h  = hs[:,:,ii]
        uh = sim.sol[sim.Iu,:,:,ii]
        vh = sim.sol[sim.Iv,:,:,ii]

        # If we're Zhang limiting, due it one way
        if sim.limiter == 'Zhang' or sim.limiter == 'zhang':
            #import Fortran_Code_Zhang
            import Fortran_Code_Zhang_2d as Fortran_Code_Zhang
            if sim.Nx > 1:
                # Step 1: Compute interpolants
                um,up = sim.axu(uh)
                vm,vp = sim.axv(vh)
                hm,hp = sim.axh(h)

                # Step 2: Compute minimum on each cell
                hmin = Fortran_Code_Zhang.zhang_minima(h,hp,hm)

                # Step 3: Compute the revised interpolants
                theta = np.min(np.dstack([np.ones(h.shape),
                                          np.abs((h-sim.min_depth)/(h - hmin))]),
                               axis=2)
                hp = theta*hp + (1-theta)*h
                hm = theta*hm + (1-theta)*h

                # Step 4: Split the flux
                Fm,Fp = Fortran_Code_Zhang.flux_split_x(h,uh,vh,hp,up,vp,hm,um,vm,sim.g)
                Fp_flux = Fp - np.roll(Fp,1,1)
                Fm_flux = np.roll(Fm,-1,1) - Fm
                the_flux[:,:,:,ii] += -0.5*(Fp_flux + Fm_flux)/sim.dx[0]
            if sim.Ny > 1:
                # Step 1: Compute interpolants
                um,up = sim.ayu(uh)
                vm,vp = sim.ayv(vh)
                hm,hp = sim.ayh(h)

                # Step 2: Compute minimum on each cell
                hmin = Fortran_Code_Zhang.zhang_minima(h,hp,hm)

                # Step 3: Compute the revised interpolants
                theta = np.min(np.dstack([np.ones(h.shape),
                                          np.abs((h-sim.min_depth)/(h - hmin))]),
                               axis=2)
                hp = theta*hp + (1-theta)*h
                hm = theta*hm + (1-theta)*h

                # Step 4: Split the flux
                Fm,Fp = Fortran_Code_Zhang.flux_split_y(h,uh,vh,hp,up,vp,hm,um,vm,sim.g)
                Fp_flux = Fp - np.roll(Fp,1,2)
                Fm_flux = np.roll(Fm,-1,2) - Fm
                the_flux[:,:,:,ii] += -0.5*(Fp_flux + Fm_flux)/sim.dx[1]

        else:
            # x fluxes
            if sim.Nx > 1:
                Fm,Fp = Fortran_Code.flux_split_x(h,uh,vh,sim.g) 

                for var in [sim.Iu,sim.Ih,sim.Iv]:
                    Fup = Fp[var,:,:]
                    Fum = Fm[var,:,:]
                    if var == sim.Iu:
                        tmp,Fxp = sim.axu(Fup)
                        Fxm,tmp = sim.axu(Fum)
                    elif var == sim.Iv:
                        tmp,Fxp = sim.axv(Fup)
                        Fxm,tmp = sim.axv(Fum)
                    elif var == sim.Ih:
                        tmp,Fxp = sim.axh(Fup)
                        Fxm,tmp = sim.axh(Fum)
                    Fp_flux = Fxp - np.roll(Fxp,1,0)
                    Fm_flux = np.roll(Fxm,-1,0) - Fxm
                    the_flux[var,:,:,ii] += -(Fp_flux + Fm_flux)/sim.dx[0]

            # y fluxes
            if sim.Ny > 1:
                Fm,Fp = Fortran_Code.flux_split_y(h,uh,vh,sim.g) 

                for var in [sim.Iu,sim.Ih,sim.Iv]:
                    Fup = Fp[var,:,:]
                    Fum = Fm[var,:,:]
                    if var == sim.Iu:
                        tmp,Fxp = sim.ayu(Fup)
                        Fxm,tmp = sim.ayu(Fum)
                    elif var == sim.Iv:
                        tmp,Fxp = sim.ayv(Fup)
                        Fxm,tmp = sim.ayv(Fum)
                    elif var == sim.Ih:
                        tmp,Fxp = sim.ayh(Fup)
                        Fxm,tmp = sim.ayh(Fum)
                    Fp_flux = Fxp - np.roll(Fxp,1,1)
                    Fm_flux = np.roll(Fxm,-1,1) - Fxm
                    the_flux[var,:,:,ii] += -(Fp_flux + Fm_flux)/sim.dx[1]

    return the_flux

# Deal with the source terms
def fv_sw_source_zhang(sim):
    
    if sim.limiter == 'zhang' or sim.limiter == 'Zhang': 
        import Fortran_Code_Zhang_2d as Fortran_Code_Zhang
    # Some storage fields
    the_source = np.zeros(sim.sol.shape)

    # Difference in layer deformations gives layer depths
    h = sim.sol[sim.Ih,:,:,:sim.Nz] - sim.sol[sim.Ih,:,:,1:]

    # Loop through each layer and compute flux

    # For layer i, the source terms look like:
    #  -g * h_i / rho_i * sum(rho_j * (d/dx)h_j) - g * h_j * (d/dx)*eta_{i+1}
    # Since the sum is built incrementally, I'll store one variable and build it as we go

    # Deal with first layer
    the_source[sim.Iu,:,:,0] +=  sim.f0*sim.sol[sim.Iv,:,:,0]
    the_source[sim.Iv,:,:,0] += -sim.f0*sim.sol[sim.Iu,:,:,0]

    hii = h[:,:,0]

    if sim.Nx > 1:
        Sxp = np.zeros((sim.Nx,sim.Ny))
        Sxm = np.zeros((sim.Nx,sim.Ny))
        
        dxm,dxp = sim.Dxh(sim.sol[sim.Ih,:,:,1],sim.dx)
        hxm,hxp = sim.axh(h[:,:,0])

        if sim.limiter == 'zhang' or sim.limiter == 'Zhang': 
            hmin = Fortran_Code_Zhang.zhang_minima(hii,hxp,hxm)
            theta = np.min(np.dstack([np.ones(hii.shape),
                                      np.abs((hii-sim.min_depth)/(hii - hmin))]),
                           axis=2)
            hxp = theta*hxp + (1-theta)*hii
            hxm = theta*hxm + (1-theta)*hii

        the_source[sim.Iu,:,:,0] += -sim.g*0.5*(dxm*hxm + dxp*hxp)
    
    if sim.Ny > 1:
        Syp = np.zeros((sim.Nx,sim.Ny))
        Sym = np.zeros((sim.Nx,sim.Ny))
        
        dym,dyp = sim.Dyh(sim.sol[sim.Ih,:,:,1],sim.dx)
        hym,hyp = sim.ayh(h[:,:,0])

        if sim.limiter == 'zhang' or sim.limiter == 'Zhang': 
            hmin = Fortran_Code_Zhang.zhang_minima(hii,hyp,hym)
            theta = np.min(np.dstack([np.ones(hii.shape),
                                      np.abs((hii-sim.min_depth)/(hii - hmin))]),
                           axis=2)
            hyp = theta*hyp + (1-theta)*hii
            hym = theta*hym + (1-theta)*hii
        
        the_source[sim.Iv,:,:,0] += -sim.g*0.5*(dym*hym + dyp*hyp)

    # Deal with rest of layers
    for ii in range(1,sim.Nz):

        hii = h[:,:,ii]

        # Deal with Coriolis terms
        the_source[sim.Iu,:,:,ii] +=  sim.f0*sim.sol[sim.Iv,:,:,ii]
        the_source[sim.Iv,:,:,ii] += -sim.f0*sim.sol[sim.Iu,:,:,ii]

        if sim.Nx > 1:
            dxm,dxp = sim.Dxh(h[:,:,ii-1],sim.dx)
            hxm,hxp = sim.axh(h[:,:,ii])
            if sim.limiter == 'zhang' or sim.limiter == 'Zhang': 
                hmin = Fortran_Code_Zhang.zhang_minima(hii,hxp,hxm)
                theta = np.min(np.dstack([np.ones(hii.shape),
                                          np.abs((hii-sim.min_depth)/(hii - hmin))]),
                               axis=2)
                hxp = theta*hxp + (1-theta)*hii
                hxm = theta*hxm + (1-theta)*hii
            Sxp += sim.rho[ii-1]*dxp
            Sxm += sim.rho[ii-1]*dxm
        
            # Deal with sum term
            the_source[sim.Iu,:,:,ii] += -(sim.g/sim.rho[ii])*0.5*(hxp*Sxp + hxm*Sxm)

            # Deal with final term
            dxm,dxp = sim.Dxh(sim.sol[sim.Ih,:,:,ii+1],sim.dx)
            the_source[sim.Iu,:,:,ii] += -sim.g*0.5*(hxp*dxp + hxm*dxm)
        
        if sim.Ny > 1:
            dym,dyp = sim.Dyh(h[:,:,ii-1],sim.dx)
            hym,hyp = sim.ayh(h[:,:,ii])
            if sim.limiter == 'zhang' or sim.limiter == 'Zhang': 
                hmin = Fortran_Code_Zhang.zhang_minima(hii,hyp,hym)
                theta = np.min(np.dstack([np.ones(hii.shape),
                                          np.abs((hii-sim.min_depth)/(hii - hmin))]),
                               axis=2)
                hyp = theta*hyp + (1-theta)*hii
                hym = theta*hym + (1-theta)*hii
            Syp += sim.rho[ii-1]*dyp
            Sym += sim.rho[ii-1]*dym
            
            # Deal with sum term
            the_source[sim.Iv,:,:,ii] += -(sim.g/sim.rho[ii])*0.5*(hyp*Syp + hym*Sym)
            
            # Deal with final term
            dym,dyp = sim.Dyh(sim.sol[sim.Ih,:,:,ii+1],sim.dx)
            the_source[sim.Iv,:,:,ii] += -sim.g*0.5*(hyp*dyp + hym*dym)

    return the_source

