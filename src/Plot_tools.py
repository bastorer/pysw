# Contains plotting commands
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

# Initialize plot objects
def initialize_plots(sim):

    fig = plt.figure()
    if sim.Nx > 1 and sim.Ny > 1:
        x = sim.x/1e3
        y = sim.y/1e3
        Qs  = []
        axs = []
        for L in range(sim.Nz):
            plt.subplot(1,sim.Nz,L+1)
            axs += [plt.gca()]
            Qs += [plt.pcolormesh(x,y,sim.sol[sim.Ih,:,:,L].T, cmap='ocean')]
            plt.clim([-1,1])
            plt.axis('tight')
            plt.gca().set_aspect('equal')

    else:
        Qs  = [[],[],[]]
        axs = []
        for var in [sim.Iu,sim.Iv,sim.Ih]:
            plt.subplot(3,1,var+1)
            axs += [plt.gca()]
            for L in range(sim.Nz):
                if sim.Nx > 1:
                    x = sim.x/1e3
                if sim.Ny > 1:
                    x = sim.y/1e3
                if var == sim.Ih:
                    l, = plt.plot(x,np.sum(sim.Hs[L:])+sim.sol[sim.Ih,:,:,L])
                    plt.ylim([0,80])
                else:
                    #if var == sim.Iu:
                    #    l, = plt.plot(x,sim.hx)
                    #if var == sim.Iv:
                    #    l, = plt.plot(x,sim.hy)
                    l, = plt.plot(x,sim.sol[var,:,:,L])
                    plt.ylim([-1.5,1.5])
                Qs[var] += [l]
                    
            if var == sim.Ih:
                # Plot topography
                plt.plot(x,sim.sol[sim.Ih,:,:,-1],'k')


    if sim.animate == 'Anim':
        plt.ion()
        plt.show()
        
    sim.Qs  = Qs
    sim.axs = axs

    if sim.animate == 'Save':
        movie_name = 'Outputs/' + sim.run_name + '.mp4'
        FFMPEG = anim.writers['ffmpeg']
        writer = FFMPEG(fps=sim.fps)
        writer.setup(fig, movie_name, sim.dpi)
        sim.movie_writer = writer

# Update plot objects
def update_plots(sim):
    for var in [sim.Iu,sim.Iv,sim.Ih]:
        for L in range(sim.Nz):
            if sim.Nx > 1 and sim.Ny > 1:
                sim.Qs[L].set_array(np.ravel(sim.sol[sim.Ih,:sim.Nx-1,:sim.Ny-1,L].T))
                sim.Qs[L].changed()
            else:
                if var == sim.Ih:
                    sim.Qs[var][L].set_ydata(np.sum(sim.Hs[L:])+sim.sol[sim.Ih,:,:,L])
                else:
                    #if var == sim.Iu:
                    #    sim.Qs[var][L].set_ydata(sim.hx)
                    #if var == sim.Iv:
                    #    sim.Qs[var][L].set_ydata(sim.hy)
                    sim.Qs[var][L].set_ydata(sim.sol[var,:,:,L])
        sim.axs[var].relim()
        sim.axs[var].autoscale_view()
        plt.draw()

    if sim.animate == 'Save':
        sim.movie_writer.grab_frame()

# Finalize
def end_movie(sim):
    sim.movie_writer.finish()

