#1D Barotropic Quasi-geostrophic model
#Tim August

#Focus on streamfunction (and Rossby wave propagation)

#Steps:
#1. Need a 1D grid 
    #1D so use an array for x values (arange)

#2. Then streamfunction perturbation 
    #Sine wave (have amp and wavenumber values)

#3. Find intial vorticty 
    #q = laplacian of streamfunction(psi)
    #So find Laplacian (finite differences from the Taylor series)

#4. Apply SOR

#5. Plots and phase speed verification


#____________________________________________________________________________________________________________________________
import numpy as np
import matplotlib.pyplot as plt

#Model constants -> keep as is
k = 2 * np.pi / lx #Wavenumber
nx = 128
lx = 1e6  #m

#Parameters -> can be modified
nt = 500
dt = 10
u = 5 #m/s
beta = 1.62e-11 #/ms
e = 0.001 #Amplitude
omega = 1.5
tol = 1e-10
t_index = 100 #This can be changed (0 to nt - 1) 
              #Also have that black line on the top plot corresponding to the time (lower plot shows snapshot of that time)


#_____________________________________________________________________________________________________________________________
#1. The grid

#Establish grid points (call nx)
#Establish lengths (call lx)
#Get grid spacing: dx 
    #dx = lx / nx
#Create an array for x

dx = lx / nx
x = np.arange(nx) * dx

#dx is grid spacing

#**Note that the quantity u*dt / dx must be below 1**
    #This ensures a proper time step
    #Keep in mind when selecting the parameters 

check = (u * dt) / dx
print("Proper time step check (<1):", check) # Below 1 and it's good

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

#The theoretical q value should simply just be -e*k**2 * sin(kx) (from laplacian operator)
#This will check it:
q_exact = -e * k**2 * np.sin(k * x)

#Then compare against the actual q
error = np.max(np.abs(q - q_exact))

print(f"Initial vorticity error = {error:.3e}")

#So, q(x) is in every grid point

#q and sine wave perturbation are now defined
#Also have a clean slate of psi values starting at 0 everywhere (from np.zeros)


#______________________________________________________________________________________________________________________________
#4. SOR
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

def SOR(q, psi, dx, omega = omega, tol = tol, max_iters = 5000): #Omega is the relaxation factor
                                                               #This is Gauss-Seidel here with omega = 1.5
                                                               #Tolerance can be played around with -> defines how small the error must be to stop iterating
    nx = len(q)
    
    for iterations in range(max_iters):
        psi_old = psi.copy()
        
        for i in range(nx):

            im = (i - 1) % nx
            ip = (i + 1) % nx
            
            psi[i] = (1 - omega) * psi[i] + (omega / 2) * (psi[im] + psi[ip] - q[i] * dx**2)
            #This comes from the second-order finite difference approximation
            #Gauss-Seidel uses new values as soon as they're available
            #The Gauss-Seidel method updates get a new value of psi[i], then the SOR moves it a little further in that direction to speed convergence
            
        #This is where tolerance is checked
        #Residual is very small
        error = np.max(np.abs(psi - psi_old))

        if error < tol:
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
    q -= np.mean(q)

    psi = SOR(q, psi, dx, omega=1.5)
    psi -= np.mean(psi)

    psi_all[n, :] = psi


#print(psi_all)

#Gonna calculate wavelength and speed here

wavelength = 2 * np.pi / k
c = u - beta / k ** 2

#print(f"Wavelength = {wavelength:.4f} m")  
#print(f"Speed = {c:.4f} m/s") 



#_____________________________________________________________________________________________________________________________
#5. Plots and phase speed verification
#Contour plot of x vs. time
#Plot of x vs. Psi
#Theoretical phase speed vs numerical phase speed verification plot
#Overlay to see results

X = x #For contour
T = np.arange(nt) * dt #For countour
                       #Gives an array for the time values to be plotted 


time_value = T[t_index] #Store the time values in this array

#Before plotting, the phase speed can be verified
#Tracking the crests of the waves can give an idea about how fast the wave is moving
#This can then be compared to the theoretical result obtained earlier in the model code here (variable: c)

#Track the crest positions
crest_positions = []

for n in range(nt):
    crest_index = np.argmax(psi_all[n, :])
    crest_positions.append(x[crest_index])

#Array for where the crest positions are at
crest_positions = np.array(crest_positions)

#Since the domain is periodic, if the wave wraps around the periodic domain, the crest tracking could become an issue
for i in range(1, len(crest_positions)):
    if crest_positions[i] - crest_positions[i-1] < -lx/2:
        crest_positions[i:] += lx
    elif crest_positions[i] - crest_positions[i-1] > lx/2:
        crest_positions[i:] -= lx

c_numerical, intercept = np.polyfit(T, crest_positions, 1)

#Results
print(f"Theoretical phase speed = {c:.4f} m/s")
print(f"Numerical phase speed = {c_numerical:.4f} m/s")


#Now can plot
plt.subplots(2, 2, layout = "compressed")

#This is for the contour fill plot
#plt.figure(figsize=(8,4))
plt.subplot(2, 2, 1)
cf = plt.contourf(X, T, psi_all, levels = np.arange(-0.0015, 0.0015, .0002), extend = "both", cmap = "rainbow") 
                                                        
#plt.colorbar(label = "Psi")

#This is the overlay of the psi profile at that time step
#Taking the black line and putting it on the contour plot, then creating the second plot for the snapshot
plt.plot(x, np.full_like(x, time_value) + psi_all[t_index, :], color = "black", lw = 1)
#Taking array for the x, t (full_like) and adding the psi plot to it
#It overlays the black line onto the contour plot; then that black line is shown below the contour plot, and you see the Rossby wave at a specified time 

plt.xlabel("x (m)")
plt.ylabel("Time (s)")
plt.xlim(0, lx)
plt.title("Psi")
#plt.show()

#This is for the contour plot with no fill
plt.subplot(2, 2, 2)
plt.contour(X, T, psi_all, levels = np.arange(-0.0015, 0.0015, .0002), cmap = "rainbow") #Gonna go with no fill here just to see another visual
#plt.colorbar(label = "Psi")

#**Note the levels = np.arange() bounds may need to be changed/modified depending on parameter input**

#Overlay the line
plt.plot(x, np.full_like(x, time_value), color = "black", lw = 1)
plt.xlabel("x (m)")
plt.ylabel("Time (s)")
plt.xlim(0, lx)
plt.title("Psi (no fill)")
#plt.show()

#This is for a plot that shows the crests for the theoretical and numerical phase speeds
plt.subplot(2, 2, 3)
x0 = crest_positions[0]

x_theory = x0 + c * T

plt.plot(T, crest_positions, color = "blue", label = "Numerical")
plt.plot(T, x_theory, color = "orange", label = "Theoretical")

plt.xlabel("Time (s)")
plt.ylabel("Crest position (m)")
plt.title("Rossby Wave Phase Propagation")
plt.legend(fontsize = 10)
#plt.show()


#This is for the plot (x vs. Psi)
plt.subplot(2, 2, 4)
plt.plot(x, psi_all[t_index, :], c = "black")
#plt.ylim(-0.002, 0.002) #These can change; just here for set conditions
plt.xlabel("x (m)")
plt.ylabel("Psi (m^2/s)")
plt.xlim(0, lx)
plt.title("Psi at a given time interval")


cbar = plt.colorbar(cf, ax = plt.gcf().axes, orientation = "horizontal")
                                        #The plt.gcf().axes just connects the color bar to all the plots
                                        #Helps give the shared color plot
                                        

cbar.set_label("Psi (m^2/s)")
plt.suptitle("Evolution of the Streamfunction (Psi)")


plt.show()

