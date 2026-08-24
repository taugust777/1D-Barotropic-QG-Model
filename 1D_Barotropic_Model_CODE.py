#1D Barotropic Quasi-geostrophic model
#Tim August

#Focus on streamfunction (and Rossby wave propagation)

#Steps:
#1. Need a 1D grid 
    #1D so use an array for x values (linspace)

#2. Then streamfunction perturbation 
    #Sine wave (have amp and wavenumber values)

#3. Find intial vorticty 
    #q = laplacian of streamfunction(psi)
    #So find Laplacian (finite differences from the Taylor series)

#4. Time step
    #Time stepping area (SOR)

#5. Plot 


#____________________________________________________________________________________________________________________________
import numpy as np
import matplotlib.pyplot as plt

#All constants 
nx = 128
lx = 1e6  #m
nt = 1000
dt = 100
u = 5 #m/s
beta = 1.62 * (10 ** -11) #/ms
e = 0.001 #Amplitude
k = 4 * np.pi / lx #Wavenumber
omega = 1.5
tol = 1e-6


#_____________________________________________________________________________________________________________________________
#1. The grid

#Establish grid points (call nx)
#Establish lengths (call lx)
#Get grid spacing: dx 
    #dx = lx / nx
#Create an array for x

dx = lx / nx
x = np.linspace(0, lx, nx)

#x from 0 to 1e6 


#dx is grid spacing


#______________________________________________________________________________________________________________________________
#2. Streamfunction perturbation

#streamfunction = psi

#Start with psi = 0 everywhere
#Introduce a sine wave
    #psi = esin(kx)
    #e is amplitude (make it small) and k is wavenumber
#This tracks instability

psi_all = np.zeros([nt, nx]) #Starts at 0

psi = e * np.sin(k * x) #This is the sine wave perturbation
#plt.plot(psi)
#plt.show()


#_______________________________________________________________________________________________________________________________
#3. Initial vorticity
#So, for 1d, q = laplacian of psi = d^2 psi / dx^2

#2nd-order finite differences (Taylor series): f double prime (x) = (f(x + delta x) - 2 * f(x) + f(x - delta x)) / delta x squared

def laplacian(f, dx):
    laplace = (np.roll(f, -1) - (2 * f) + np.roll(f, 1)) / (dx ** 2) #This is the second-order finite difference approx (Taylor series)
    return laplace
    #The np.roll function is basically just taking the old function, subtracting double the function, and then adding the new function
    #so like np.roll(f, -1) shifts all to the left
            #np.roll(f, 1) shifts all to the right
    #(fi + 1) − 2 * fi + (fi − 1) -> Its like wrapping around 

q = laplacian(psi, dx) #This is the initial q
#print(q)

#plt.plot(q)
#plt.show()

#So, q(x) is in every grid point

#q and sine wave perturbation are now defined
#Also have a clean slate of psi values starting at 0 everywhere (from np.zeros)


#______________________________________________________________________________________________________________________________
#4. Time step
#Use dq/dt + u * dq/dx + beta * d psi/dx = 0
#Then dq/dt = -u * dq/dx - beta * d psi/dx

#Set u and beta, then just calculate dq/dx and d psi/dx while in range

#u = 5 m/s -> mid-lat troposphere assumption
#beta = 1.62 * (10 ** -11) /ms -> Working aorund 45 degrees latitude here

#These values, u and beta, can change
#Found some typical values for these
#u varies depending on level or location
    #mid-lat, troposphere about 5-15 m/s
    #Jet stream about 20-50 m/s
    #Stratospheric jet about 10-40 m/s
#Beta depends on phi (latitude) -> beta = df/dy = 1/a * d/dphi * 2 * omega * sin(phi) = 2 * omega * cos(phi) / a
    #phi of 0 degrees, beta ~ 2.29 * 10^-11 /ms
    #phi of 30 degrees, beta ~ 1.98 * 10^-11 /ms
    #phi of 45 degrees, beta ~ 1.62 * 10^-11 /ms
    #phi of 60 degrees, beta ~ 1.14 * 10^-11 /ms

#Need to establish a time frame
#Step forward by dt
#Number of time frames will be nt


#total time = nt * dt

#This is the successive over-relaxation method (SOR):

#The method summary:
    #Have a current guess
    #Compute residual (or how far you are from satisfying the equation)
    #The residual is used to correct the guess (scaled by an omega term)
    #omega > 1 takes bigger corrective steps -> over relaxation
    #Sweep across all iterative values
    #Repeat until residuals are below tolerance

#Matrix form SOR: (D + omega*L)x (k+1) = omega*b - [omega*U + (omega-1)D]x (k)
def SOR(q, dx, omega = 1.5, tol = 10**(-6), max_iters = 5000): #Omega is the relaxation factor
                                                               #This is Gauss-Seidel here with omega = 1.5
                                                               #Tolerance can be played around with -> defines how small the error must be to stop iterating
    nx = len(q)
    for iterations in range(max_iters):
        psi_old = psi.copy()
        
        for i in range(1, nx - 1):
            psi[i] = (1 - omega) * psi[i] + (omega / 2) * (psi[i-1] + psi[i+1] + q[i] * dx**2)
            #This comes from the second-order finite difference approximation
            #Gauss-Seidel uses new values as soon as they're available
            #The Gauss-Seidel method updates get a new value of psi[i], then the SOR moves it a little further in that direction to speed convergence
            #The weighted average is taken, which is what introduces the omega relaxation term

        #This is where tolerance is checked
        #Residual is very small
        if np.max(np.abs(psi - psi_old)) < tol:
            break
    return psi

#print(psi)

#plt.plot(psi)
#plt.show()

#Now use the for loop
for n in range(nt):
    dq_dx = (np.roll(q, -1) - np.roll(q, 1)) / (2 * dx)
    dpsi_dx = (np.roll(psi, -1) - np.roll(psi, 1)) / (2 * dx)
    #This is from the taylor series finite differences method for first derivatives
    #f prime (x) = (f(x + delta x) - f(x - delta x)) / 2 * delta x

    #Then get the time derivative dq/dt
    dq_dt = -u * dq_dx - beta * dpsi_dx

    #Then step it forward
    q += dt * dq_dt

    psi = SOR(q, dx, omega = 1.5)
    psi_all[n, :] = psi

    #Trying to get psi (to plot it)
    #Need to update psi
    #That comes from q
    #And q is the Laplacian of psi
    #Use successive over-relaxation - Gauss-Seidel method
    #Gonna just start with omega = 1 then work up from there


#print(psi_all)

#Gonna calculate wavelength and speed here

wavelength = 2 * np.pi / k
c = u - beta / k ** 2

print(f"Wavelength = {wavelength:.4f} m")  #Based on the assigned k, it should be about half of lx (1 * 10^6)
print(f"Speed = {c:.4f} m/s") #Should be close to u value at the start



#_____________________________________________________________________________________________________________________________
#5. Plots
#Contour plot of x vs. time
#Plot of x vs. Psi
#Overlay to see results

X = x #For contour
T = np.arange(nt) * dt #For countour
                       #Gives an array for the time values to be plotted 

t_index = 200 #This can be changed (0 - 999) -> multiplied by 10^2 on y axis 
              #Also have that black line on the top plot corresponding to the time (lower plot shows snapshot of that time)
time_value = T[t_index] #Store the time values in this array


#This is for the contour plot with 25 levels
#plt.figure(figsize=(8,4))
plt.subplot(2, 2, 1)
cf = plt.contourf(X, T, psi_all, levels = 25, cmap = "rainbow") 
                                                        
#plt.colorbar(label = "Psi")

#This is the overlay of the psi profile at that time step
#Taking the black line and putting it on the contour plot, then creating the second plot for the snapshot
plt.plot(x, np.full_like(x, time_value) + psi_all[t_index, :], color = "black", lw = 1)
#Taking array for the x, t (full_like) and adding the psi plot to it
#It overlays the black line onto the contour plot; then that black line is shown below the contour plot, and you see the Rossby wave at a specified time 

plt.xlabel("x (m)", labelpad = -1)
plt.ylabel("Time (s)")
plt.title("Psi (25 levels)", y = 0.98)
#plt.show()

#This is for the contour plot with no fill
plt.subplot(2, 2, 2)
plt.contour(X, T, psi_all, cmap = "rainbow") #Gonna go with no fill here just to see another visual
#plt.colorbar(label = "Psi")

#Overlay the line
plt.plot(x, np.full_like(x, time_value), color = "black", lw = 1)
plt.xlabel("x (m)", labelpad = -1)
plt.ylabel("Time (s)")
plt.title("Psi (no fill)", y = 0.98)
#plt.show()

#This is for the contour plot with 50 levels
plt.subplot(2, 2, 3)
plt.contourf(X, T, psi_all, levels = 50, cmap = "rainbow") #Gonna go with 50 levels and fill here just to see another visual
#plt.colorbar(label = "Psi")

#Overlay the line
plt.plot(x, np.full_like(x, time_value) + psi_all[t_index, :], color = "black", lw = 1)
plt.xlabel("x (m)")
plt.ylabel("Time (s)")
plt.title("Psi (50 levels)", y = 0.98)
#plt.show()

#This is for the plot (x vs. Psi)
plt.subplot(2, 2, 4)
plt.plot(x, psi_all[t_index, :], c = "black")
plt.ylim(-0.0025, 0.0025) #These can change; just here for set conditions
plt.xlabel("x (m)")
plt.ylabel("Psi (m^2/s)")
plt.title("Psi at a given time interval", y = 0.98)


plt.subplots_adjust(hspace = 0.6, wspace = 0.3) #This is for spacing between plots

#Adjusting the layout and adding shared colorbar here
plt.tight_layout(rect = [0, 0.05, 1, 1]) #This rect function is adding space for that common/shared color bar
cbar = plt.colorbar(cf, ax = plt.gcf().axes, orientation = "horizontal", fraction=0.05, pad=0.17)
                                        #The plt.gcf().axes just connects the color bar to all the plots
                                        #Helps give the shared color plot
                                        #Fraction is the width of the shared color bar relative to my plots
                                        #Pad moves the location of my shared color bar
cbar.set_label("Psi (m^2/s)")
plt.suptitle("Evolution of the Streamfunction (Psi)", y = 1.05)
#plt.legend(print(f"Wavelength = {wavelength:.4f} m"), print(f"Speed = {c:.4f} m/s"))
plt.figure(figsize = (10, 14))
plt.show()

