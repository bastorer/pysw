import numpy as np
import matplotlib.pyplot as plt
import FV_plot_tools as Plot_tools
import FV_Diagnose as Diagnose
import os
import shutil

class Sadourny:

    # First-level initialization
    def __init__(self):
        # Default values
        self.nfluxes = 0        # Length of flux history
        self.Iu = 0             # u index
        self.Iv = 1             # v index
        self.Ih = 2             # h index
        self.fluxes = []        # flux history
        self.g    = 9.81        # gravity
        self.f0   = 1e-4        # Coriolis
        self.cfl  = 0.01        # default CFL
        self.time = 0           # initial time
        self.min_depth = 0.1    # minimum fluid depth
        self.vanishing = False
        self.fps  = 15          # movie fps
        self.dpi  = 150         # movie dpi
        self.run_name = 'test'  # movie name
        self.min_dt = 1e-10     # minimum timestep
        self.BCx = 'periodic'   # x boundary condition
        self.BCy = 'periodic'   # y boundary condition
        self.ylims = [[],[],[]]

    # Initialize everything else
    def initialize(self):

        # Initialize grids (using Arakawa C grid)
        dxs = [1,1]
        if self.Nx > 1:
            dx = self.Lx/self.Nx
            self.x = np.arange(dx/2,self.Lx,dx)
            dxs[0] = dx
            print('dx = ',dx)
        if self.Ny > 1:
            dy = self.Ly/self.Ny
            self.y = np.arange(dy/2,self.Ly,dy)
            dxs[1] = dy
            print('dy = ', dy)
        self.dx = dxs

        # Initialize differentiation and averaging operators
        self.x_derivs(self)
        self.y_derivs(self)

        # Compute reduced gravities
        if self.Nz > 1:
            tmp = np.ravel(self.rho.copy())
            tmp[1:] = tmp[1:] - tmp[0:-1]
            tmp = np.tile(tmp,(1,1)).reshape((1,len(self.rho)))
            self.gs = self.g*tmp/tmp[0,0]
        else:
            self.gs = np.array([[self.g]])

        # Initial conditions and topography
        self.sol = np.zeros((3,self.Nx,self.Ny,self.Nz+1)) # Extra layer for topography
        self.topo_func(self)
        self.IC_func(self)

        # If we're going to be plotting, then initialize the plots
        if self.animate != 'None':
            Plot_tools.initialize_plots(self)
            self.next_plot_time = self.ptime

        # If we're going to be diagnosing, initialize those
        Diagnose.initialize_diagnostics(self)
    
        # If we're saving, initialize those too
        if self.output:
            self.out_counter = 0
            self.initialize_saving()
            self.next_save_time = self.otime

    # Compute the current flux
    def flux(self):
        return self.flux_function(self)

    # Compute the current sources
    def source(self):
        return self.source_function(self)

    # Adjust time-step when necessaru
    def adjust_dt(self):
        
        t = self.time + self.dt

        nt = np.Inf # next time for doing stuff
        if self.animate != 'None':
            nt = min([self.next_plot_time, nt])
        if self.diagnose:
            nt = min([self.next_diag_time, nt])
        if self.output:
            nt = min([self.next_save_time, nt])

        if nt < t:
            self.dt = t - nt

        t = self.time + self.dt

        do_plot = False
        if self.animate == 'Anim' or self.animate == 'Save':
            if t + self.min_dt >= self.next_plot_time:
                do_plot = True

        do_diag = False
        if self.diagnose:
            if t + self.min_dt >= self.next_diag_time:
                do_diag = True

        do_save = False
        if self.output:
            if t + self.min_dt >= self.next_save_time:
                do_save = True


        return do_plot, do_diag, do_save

    # Advance the simulation one time-step.
    def step(self):
        self.compute_dt()

        # Check if we need to adjust the time-step
        # to match an output time
        do_plot, do_diag, do_save = self.adjust_dt()

        self.time_stepper(self)
        self.time += self.dt
       
        if do_plot:
            Plot_tools.update_plots(self)
            self.next_plot_time += self.ptime

        if do_diag:
            Diagnose.update(self)

        if do_save:
            self.save_state()
            self.next_save_time += self.otime

        if do_plot or self.time == self.dt:

            h0 = self.sol[self.Ih,:,:,0] - self.sol[self.Ih,:,:,1]
            minh = np.min(np.ravel(h0))
            maxh = np.max(np.ravel(h0))
            maxu = np.max(np.ravel(self.sol[self.Iu,:,:,0]/h0))
            minu = np.min(np.ravel(self.sol[self.Iu,:,:,0]/h0))
            maxv = np.max(np.ravel(self.sol[self.Iv,:,:,0]/h0))
            minv = np.min(np.ravel(self.sol[self.Iv,:,:,0]/h0))
            mass = Diagnose.compute_mass(self)
            enrg = Diagnose.compute_PE(self) + Diagnose.compute_KE(self)

            if self.Nz == 1:
                pstr  = 't = {0: <5.4g}s'.format(self.time)
                pstr += ' '*(13-len(pstr))
                try:
                    pstr += ', dt = {0:0<7.1e}'.format(np.mean(self.dts))
                except:
                    pstr += ', dt = {0:0<7.1e}'.format(self.dt)
                L = len(pstr) - 11
                pstr += ', min(u,v,h) = ({0: < 8.4e},{1: < 8.4e},{2: < 8.4e})'.format(minu,minv,minh)
                pstr += ', del_mass = {0: .2g}'.format(mass/self.Ms[0]-1)
                pstr += '\n'
                tmp = '  = {0:.3%}'.format(self.time/self.end_time)
                pstr += tmp
                pstr += ' '*(L + 13 - len(tmp))            
                pstr += 'max(u,v,h) = ({0: < 8.4e},{1: < 8.4e},{2: < 8.4e})'.format(maxu,maxv,maxh)
                pstr += ', del_enrg = {0: .2g}'.format(enrg/(self.KEs[0]+self.PEs[0])-1)
            else:
                pstr  = 't = {0:.4g}s'.format(self.time)
                try:
                    pstr += ', dt = {0:0<7.1e}'.format(np.mean(self.dts))
                except:
                    pstr += ', dt = {0:0<7.1e}'.format(self.dt)
                pstr += ', del_mass = {0:+.2g}'.format(mass/self.Ms[0]-1)
                pstr += ', del_enrg = {0:+.2g}'.format(enrg/(self.KEs[0]+self.PEs[0])-1)
                pstr += '\n'
                pstr += '  = {0:.3%}'.format(self.time/self.end_time)
            print('\n{0:s}: {1:d}'.format(self.run_name, int(self.time/self.ptime)))
            print(pstr)

    # Step until end-time achieved.
    def run(self):

        while self.time < self.end_time:
            self.step()

        if self.diagnose:
            Diagnose.plot(self)
        
        if self.animate == 'Anim':
            plt.ioff()
            plt.show()

        if self.animate == 'Save':
            Plot_tools.end_movie(self)

        if self.diagnose:
            Diagnose.save(self)

    # Compute time-step using CFL condition.
    def compute_dt(self):
        c = np.sqrt(self.g*np.sum(self.Hs))
        eps = 1.e-8 
        hs = self.sol[self.Ih,:,:,:self.Nz] - self.sol[self.Ih,:,:,1:]

        if self.vanishing:
            tmp = self.min_depth**self.np/(hs**(1-self.np))
            c = np.max(np.sqrt(c**2 + self.g*tmp).ravel())

        if self.Nx > 1:
            u = self.sol[self.Iu,:,:,:self.Nz]/(eps + hs)
            max_u = np.max(np.abs(u.ravel()))
            dt_x = self.dx[0]/(max_u+2*c)
        else:
            dt_x = np.Inf

        if self.Ny > 1:
            u = self.sol[self.Iv,:,:,:self.Nz]/(eps + hs)
            max_v = np.max(np.abs(u.ravel()))
            dt_y = self.dx[1]/(max_v+2*c)
        else:
            dt_y = np.Inf

        self.dt = max([self.cfl*min([dt_x,dt_y]),self.min_dt])

    # Initialize the saving
    def initialize_saving(self):
        
        path = 'Outputs/{0:s}'.format(self.run_name)
        # If directory already exists, delete it.
        if os.path.isdir(path):
           print('Output directory {0:s} already exists. '.format(path) + \
                 'Warning, deleting everything in the directory.') 
           shutil.rmtree(path)

        # Make directory.
        os.mkdir(path)

        # Initialize saving stuff
        self.save_state()
        self.save_info()
        self.save_grid()

    # Save information
    def save_info(self):
        fname = 'Outputs/{0:s}/info'.format(self.run_name)

        fp = open(fname, 'w')

        fp.write('Hs         = {0:s}\n'.format(self.array2str(self.Hs)))
        fp.write('rho        = {0:s}\n'.format(self.array2str(self.rho)))
        fp.write('g0         = {0:g}\n'.format(self.g))
        fp.write('Nx         = {0:d}\n'.format(self.Nx))
        fp.write('Ny         = {0:d}\n'.format(self.Ny))
        fp.write('Nz         = {0:d}\n'.format(self.Nz))
        fp.write('min_dt     = {0:g}\n'.format(self.min_dt))
        fp.write('min_depth  = {0:g}\n'.format(self.min_depth))

        fp.close()

    # Save grid 
    def save_grid(self):
        fname = 'Outputs/{0:s}/grid'.format(self.run_name) 
        if self.Nx > 1 and self.Ny > 1:
            np.savez_compressed(fname, x = self.x, y = self.y, xs = self.xs, ys = self.ys)
        elif self.Nx > 1:
            np.savez_compressed(fname, x = self.x, xs = self.xs)
        elif self.Ny > 1:
            np.savez_compressed(fname, y = self.y, ys = self.ys)

    # Save current state
    def save_state(self):
        np.savez_compressed('Outputs/{0:s}/{1:04d}'.format(self.run_name,self.out_counter),
                    sol = self.sol, t = self.time)

        self.out_counter += 1

    # Helper to write arrays
    def array2str(fp,arr):
        s = '[' + reduce(lambda s1,s2: s1+', '+s2, map(str,arr)) + ']'
        return s

