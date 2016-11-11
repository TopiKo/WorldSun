import matplotlib.pyplot as plt
import matplotlib.animation as animation
import numpy as np
import matplotlib as mpl
from sunShines import functions, sunShine




class AnimatedSunShine(object):
    """An animated scatter plot using matplotlib.animations.FuncAnimation."""
    def __init__(self, n = 120, ndays = 365, lat_d = 62, lon_d = 0):

        self.n = n
        self.ndays = ndays
        self.get_data(lat_d, lon_d)
        interval = int(1000*30./len(self.data))

        #self.data = data
        #self.times = times
        #self.riseAndSet_times = riseAndSet
        #self.fracs = fracs
        # Setup the figure and axes...
        self.fig, self.ax = plt.subplots(figsize = (10,10))
        #self.time_text = self.ax.text(-1, 1., 'test') #, transform=self.ax.transAxes)

        # Then setup FuncAnimation.
        self.ani = animation.FuncAnimation(self.fig, self.update, frames = len(self.data),
                                           interval=interval, init_func=self.setup_plot,
                                           blit=True, repeat = False)


    def get_data(self, lat_d, lon_d):
        #n = 120
        #ndays = 365
        n, ndays = self.n, self.ndays
        funcs = functions(lat_d = lat_d, lon_d = lon_d)
        sun_riseAndSet, energies = [], []
        dayCount = 3
        #data,fracs,times = sunShine(funcs, 0, np.linspace(0,24, n))[:3]
        for i in range(int(ndays/dayCount)):
            i *= dayCount
            read_data = sunShine(funcs, i, np.linspace(0,24, n))[:4]
            if read_data[0] != None:
                day = int(read_data[2][0]/24)
                times_use = read_data[2]/24 - int(read_data[2][0]/24)

                if 'data' in locals():
                    data = np.concatenate((data, read_data[0]))
                    fracs = np.concatenate((fracs,read_data[1]))
                    times = np.concatenate((times,read_data[2]))
                    energy = read_data[3]
                if 'data' not in locals() and read_data[0] != None:
                    times_use = read_data[2]/24 - day
                    data,fracs,times,energy = read_data
                energies.append([day, energy])
                sun_riseAndSet.append([day, np.min(times_use)*24, np.max(times_use)*24])

        self.data = data
        self.fracs = fracs
        self.times = times
        self.riseAndSet_times = np.array(sun_riseAndSet)
        self.energies = np.array(energies)

    def setup_plot(self):
        """Initial drawing of the scatter plot."""
        # PLot init lines for projection
        phi = np.linspace(0, 2*np.pi, 100)
        thetas = np.linspace(np.pi/2, 0, 6, endpoint = True)

        self.n = mpl.colors.Normalize(vmin = min(self.fracs), vmax = max(self.fracs))
        self.m = mpl.cm.ScalarMappable(norm=self.n, cmap=mpl.cm.afmhot)

        for theta in thetas:
            xs,ys,zs = np.sin(theta)*np.cos(phi), np.sin(theta)*np.sin(phi), -np.cos(theta)
            xprojs, yprojs = xs/(1-zs), ys/(1-zs)
            self.ax.plot(xprojs, yprojs, color = 'black')
        for phi in np.linspace(0, 2*np.pi, 8, endpoint = False):
            xs,ys,zs = np.sin(thetas)*np.sin(-phi), np.sin(thetas)*np.cos(-phi), -np.cos(thetas)
            xprojs, yprojs = xs/(1-zs), ys/(1-zs)
            self.ax.plot(xprojs, yprojs, color = 'black')
            self.ax.text(1.1*np.sin(phi) - .05, 1.1*np.cos(-phi) - .02,
                '%.0f' %(phi/np.pi*180), fontsize = 15)

        self.ax.axes.get_xaxis().set_visible(False)
        self.ax.axes.get_yaxis().set_visible(False)
        self.ax.axis('equal')
        self.time_text = self.ax.text(-1., .82, 'test', fontsize = 15)

        circle = plt.Circle((0, 0), 1.0, color=(40./255, 150./255, 230./255))
        self.ax.add_artist(circle)
        x,y = self.data[0]

        self.scat = self.ax.scatter(x, y, animated=True, zorder = 10000)
        self.ax.axis([-1, 1, -1, 1])
        plt.axis('off')
        # For FuncAnimation's sake, we need to return the artist we'll be using
        # Note that it expects a sequence of artists, thus the trailing comma.
        return self.scat, self.time_text

    def update(self, i):
        """Update the scatter plot."""
        data = self.data
        num = i
        if 10 < num:
            xs = data[num - 10:num,0]
            ys = data[num - 10:num,1]
            fracs = np.array(self.fracs[num-10:num])
        else:
            xs = data[:num,0]
            ys = data[:num,1]
            fracs = np.array(self.fracs[:num])

        plotPos = np.zeros((len(xs),2))
        plotPos[:,0] = xs
        plotPos[:,1] = ys

        time = self.times[num]/24
        day = int(time)
        #k = int(day/dayCount)
        try:
            k = np.where(self.riseAndSet_times[:,0] == day)[0][0]
            self.last_k = k
        except IndexError as ie:
            k = self.last_k
        hour = int(time*24 - day*24)
        #print(int(self.riseAndSet_times[k][0]))
        sun_rise_h = int(self.riseAndSet_times[k][1])
        sun_rise_min = int((self.riseAndSet_times[k][1] - sun_rise_h)*60)
        sun_set_h = int(self.riseAndSet_times[k][2])
        sun_set_min = int((self.riseAndSet_times[k][2] - sun_set_h)*60)

        self.time_text.set_text('day %03d \n' %int(time) + 'rise: ' +
        '%02d:%02d \n' %(sun_rise_h, sun_rise_min) + 'set: %02d:%02d' %(sun_set_h, sun_set_min))

        if len(xs) >2:
            self.scat.set_offsets(plotPos)
            self.scat.set_facecolor(self.m.to_rgba(fracs))
            self.scat.set_edgecolor(self.m.to_rgba(fracs))

            self.scat.set_clim(vmin=min(fracs), vmax=max(fracs))

            self.scat._sizes = 1200 * fracs**2
        return self.scat, self.time_text

        # Set colors..
        #self.scat.set_array(data[3])
        # We need to return the updated artist for FuncAnimation to draw..
        # Note that it expects a sequence of artists, thus the trailing comma.
        #return self.scat,

    def show(self):
        plt.show()

    def write(self):
        Writer = animation.writers['ffmpeg']
        writer = Writer(fps=35, metadata=dict(artist='Me'), bitrate=1800)
        self.ani.save('im.mp4', writer=writer)

    def plot_energies(self, lats):
        plt.close('all')

        for lat_d in lats:
            self.get_data(lat_d, 0)
            plt.plot(self.energies[:,0], self.energies[:,1],
                    label = 'lat deg = %.0f' %lat_d)
        plt.legend(loc = 1)
        plt.xlabel('day of year')
        plt.ylabel('radiation energy on ground in a day')
        plt.savefig('radiation_energy.png')
        plt.show()

if __name__ == '__main__':
    a = AnimatedSunShine(lat_d = 61.6)
    #a.show()
    a.plot_energies([0, 20, 40, 60])
    #a.write()
