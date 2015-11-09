# Contains plotting commands
import matplotlib.pyplot as plt
import matplotlib.animation as anim
import numpy as np

# Helper for formatting time strings
def smart_time(t):
    tstr = 't = '

    if t < 2*60.:
        tstr += '{0:.4f} sec'.format(t)
    elif t < 90.*60:
        tstr += '{0:.4f} min'.format(t/60)
    elif t < 48.*60*60:
        tstr += '{0:.4f} hrs'.format(t/(60*60))
    else:
        tstr += '{0:.4f} day'.format(t/(24*60*60))

    return tstr

# Initialize plot objects
def initialize_plots(sim):

    fig = plt.figure()
    sim.fig = fig
    fig.suptitle('t = 0')

    if sim.Nx > 1 and sim.Ny > 1:
        x = sim.x/1e3
        y = sim.y/1e3
        Qs  = []
        axs = []
        for L in range(sim.Nz):
            plt.subplot(1,sim.Nz,L+1)
            axs += [plt.gca()]
            Qs += [plt.pcolormesh(x,y,sim.sol[sim.Ih,:,:,L].T, cmap='ocean')]
            plt.colorbar()
            plt.contour(x,y,sim.sol[sim.Ih,:,:,-1].T)
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
                   
                l, = plt.plot(x,sim.sol[var,:,:,L].ravel())
                if len(sim.ylims[var]) == 2:
                    plt.ylim(sim.ylims[var])
                Qs[var] += [l]
                    
            if var == sim.Ih:
                # Plot topography
                plt.plot(x,sim.sol[sim.Ih,:,:,-1].ravel(),'k')


    if sim.animate == 'Anim':
        plt.ion()
        plt.show()
        
    sim.Qs  = Qs
    sim.axs = axs

    if sim.animate == 'Save':
        pass
        #movie_name = 'Outputs/' + sim.run_name + '.mp4'
        #FFMPEG = anim.writers['ffmpeg']
        #writer = FFMPEG(fps=sim.fps)
        #writer.setup(fig, movie_name, sim.dpi)
        #sim.movie_writer = writer
        #sim.movie_writer.grab_frame()

# Update plot objects
def update_plots(sim):

    sim.fig.suptitle(smart_time(sim.time))

    if sim.Nx > 1 and sim.Ny > 1:
        for L in range(sim.Nz):
            sim.Qs[L].set_array(np.ravel(sim.sol[sim.Ih,:sim.Nx-1,:sim.Ny-1,L].T))
            sim.Qs[L].changed()
    else:
        for var in [sim.Iu,sim.Iv,sim.Ih]:
            for L in range(sim.Nz):
                sim.Qs[var][L].set_ydata(sim.sol[var,:,:,L])
            sim.axs[var].relim()
            sim.axs[var].autoscale_view()
    plt.draw()

    if sim.animate == 'Save':
        sim.fig.savefig('Frames/{0:04d}.png'.format(sim.frame_count))
        sim.frame_count += 1

# Finalize
def end_movie(sim):
    pass
    #sim.movie_writer.finish()

