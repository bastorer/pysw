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
            Qs += [plt.pcolormesh(x,y,sim.soln.h[:,:,L].T, cmap='viridis')]
            plt.colorbar()
            try:
                plt.contour(x,y,sim.soln.h[:,:,-1].T)
            except:
                pass
            plt.axis('tight')
            plt.gca().set_aspect('equal')

    else:
        Qs  = [[],[],[]]
        axs = []

        # Plot u
        plt.subplot(3,1,1)
        axs += [plt.gca()]
        for L in range(sim.Nz):
            if sim.Nx > 1:
                x = sim.x/1e3
            if sim.Ny > 1:
                x = sim.y/1e3
            l, = plt.plot(x,sim.soln.u[:,:,L].ravel())
            if len(sim.ylims[0]) == 2:
                plt.ylim(sim.ylims[0])
            if sim.method == 'Spectral':
                plt.ylabel('u')
            else:
                plt.ylabel('uh')
            Qs[0] += [l]

        # Plot v
        plt.subplot(3,1,2)
        axs += [plt.gca()]
        for L in range(sim.Nz):
            if sim.Nx > 1:
                x = sim.x/1e3
            if sim.Ny > 1:
                x = sim.y/1e3
            l, = plt.plot(x,sim.soln.v[:,:,L].ravel())
            if len(sim.ylims[1]) == 2:
                plt.ylim(sim.ylims[1])
            if sim.method == 'Spectral':
                plt.ylabel('v')
            else:
                plt.ylabel('vh')
            Qs[1] += [l]

        # Plot h
        plt.subplot(3,1,3)
        axs += [plt.gca()]
        for L in range(sim.Nz):
            if sim.Nx > 1:
                x = sim.x/1e3
            if sim.Ny > 1:
                x = sim.y/1e3
            l, = plt.plot(x,sim.soln.h[:,:,L].ravel())
            if len(sim.ylims[2]) == 2:
                plt.ylim(sim.ylims[2])
            plt.ylabel('h')
            Qs[2] += [l]

        plt.plot(x,sim.soln.h[:,:,-1].ravel(),'k')


    if sim.animate == 'Anim':
        plt.ion()
        plt.pause(0.01)
        plt.draw()
        
    sim.Qs  = Qs
    sim.axs = axs


# Update plot objects
def update_plots(sim):

    sim.fig.suptitle(smart_time(sim.time))

    if sim.Nx > 1 and sim.Ny > 1:
        for L in range(sim.Nz):
            sim.Qs[L].set_array(np.ravel(sim.soln.h[:sim.Nx-1,:sim.Ny-1,L].T))
            sim.Qs[L].changed()
    else:
        # Update u
        for L in range(sim.Nz):
            sim.Qs[0][L].set_ydata(sim.soln.u[:,:,L])
        sim.axs[0].relim()
        sim.axs[0].autoscale_view()

        # Update v
        for L in range(sim.Nz):
            sim.Qs[1][L].set_ydata(sim.soln.v[:,:,L])
        sim.axs[1].relim()
        sim.axs[1].autoscale_view()

        # Update h
        for L in range(sim.Nz):
            sim.Qs[2][L].set_ydata(sim.soln.h[:,:,L])
        sim.axs[2].relim()
        sim.axs[2].autoscale_view()

    if sim.animate == 'Anim':
        plt.pause(0.01)
    plt.draw()

    if sim.animate == 'Save':
        sim.fig.savefig('Frames/{0:05d}.png'.format(sim.frame_count))
        sim.frame_count += 1

# Finalize
def end_movie(sim):
    pass
