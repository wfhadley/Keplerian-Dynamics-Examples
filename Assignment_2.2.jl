### A Pluto.jl notebook ###
# v0.19.37

using Markdown
using InteractiveUtils

# ╔═╡ fd9d5608-8d6f-4904-af06-c93d5b9a05bd
using PlutoUI, LinearAlgebra, Plots

# ╔═╡ e052fe10-c756-4128-b5e4-7c19d2d389a6
md"""
# 2.2
"""

# ╔═╡ 04db8de6-4267-4df9-b38a-217a46cce879
begin
	# Parameters
	
	μₛ = 1.3271244e11 # Sun gravitational parameter [km³/s²]
	μₑ = 3.986e5 # Earth gravitational parameter [km³/s²]
	Rₑ = 6378.14 # Earth radius [km]
	J₂ = 1082.63e-6 # Earth J₂ Perturbation
end;

# ╔═╡ 44993ac9-2ed3-4f05-8fdf-245cee136561
# Initial State Vector
sᵢ = [−3174.89878864632, −992.843201535607, −6922.86051061034, 5.29933974840672, 3.41395143650654, −2.9199470061499];

# ╔═╡ c38eae7c-b9f5-4591-83ce-4e1f6f04c8ed
md"""
### 2.2.1
A function is created which converts the state vector $(r, v)$ of a spacecraft from Earth Centered Inertial (ECI) reference position and velocity to Classical Orbital Elements (COE). This will allow the conversion of the initial state vector:

```
sᵢ = [x, y, z, u, v, w] = [−3174.89878864632, −992.843201535607, −6922.86051061034, 5.29933974840672, 3.41395143650654, −2.9199470061499]
```

from ECI reference coordinates to COE. The process for this conversion is as follows:
- Calculate the magnitudes of $\vec{r}$ and $\vec{v}$ ($r$ and $v$)
- Calculate the radial velocity $v_r = \frac{\vec{r} \cdot \vec{v}}{r}$. Note that if $v_r > 0$ then the positional radius is increasing, which implies that the spacecraft is moving towards apocentre. Correspondingly, if $v_r < 0$, then the spacecraft is moving towards pericentre.
- Calculate the specific angular momentum as:
$\vec{h} = \vec{r} \times \vec{v}$
- Take the magnitude of $\vec{h}$:
$\|\vec{h}\| = \sqrt{\vec{h} \cdot \vec{h}}$
- Compute the semi-major axis, $a$, from the Vis-Viva equation:
$\frac{1}{2} v^2 - \frac{\mu}{r} = -\frac{\mu}{2 a}$
$\implies a = \mu \frac{1}{2(\frac{\mu}{r} - \frac{1}{2} v^2)}$
- Compute the eccentricity vector using the definition for 2-body Keplerian dynamics:
$\vec{e} = \frac{\vec{\dot{r}} \times \vec{h}}{\mu} - \frac{\vec{r}}{r}$
This vector has constant direction and magnitude, and lies in the orbital plane ($\vec{e} \cdot \vec{h} = 0$).
- Compute the inclination as:
$i = \arccos(\frac{h_z}{h})$
This formula comes from the dot product - if it is considered that $\vec{h}$ is defined in the ECI reference frame, and $i$ is the angle between $\vec{h}$ (direction normal to orbital plane) and $\hat{k}$ (direction normal to to the $XY$ plane in the ECI frame), then it follows that the third component of $\vec{h}$ is in the $\hat{k}$ direction, and the dot product formula may be reformulated to find the angle $i$:

$\|\vec{h}\|\|\hat{k}\|cos(i) = \vec{h} \cdot \hat{k} = h_z$
$\implies i = \arccos(\frac{h_z}{h})$

Note that $i$ must take values on the interval [0°, 180°], which matches the range of the inverse cosine function.
- Calculate the line of nodes using the definition:
$\vec{n} = \hat{k} \times \vec{h}$
- Calculate the magnitude of $\vec{n}$:

$\|\vec{n}\| = \sqrt{\vec{n} \vec{n}}$
- Using a similar method to that used for the inclination, the longitude of ascending node is computed as:
$\Omega = \arccos{\frac{n_x}{n}}$
This equation gives a value on the interval [0°, 180°], however, Ω takes values on the interval [0°, 360°]. As such, a method is constructed to place Ω in the correct quadrant of its plane. To place Ω in the correct quadrant, the $\hat{j}$ component of $\vec{n}$ is observed. If this value is positive, then Ω lies in quadrant 1 or 2 of the $\hat{i}, \hat{j}$ plane. If this value is negative, then the quadrant is either 3 or 4. In summary:

>$\begin{align}
>\Omega&=\arccos{\frac{n_x}{n}}, (n_y \ge 0),\\
>\Omega&= 360° - \arccos{\frac{n_x}{n}}, (n_y < 0)
>\end{align}$

- Calculate the argument of perigee using a similar method to the calculation of $i$ and $\Omega$:

$\omega  = \arccos{\frac{\vec{n} \cdot \vec{e}}{n e}}$

To place $\omega$ in the correct quadrant, observe that perigee lies above the equatorial plane if $\vec{e}$ points up (in positive z direction) and that perigee lies below the plane if $\vec{e}$ points down (negative z direction). In summary:

>$\begin{align}
>\omega&=\arccos{\frac{\vec{n} \cdot \vec{e}}{n e}}, (e_z \ge 0),\\
>\omega&= 360° - \arccos{\frac{\vec{n} \cdot \vec{e}}{n e}}, (e_z < 0)
>\end{align}$

- Compute the true anomaly by the same dot product reformulation method as has been used above:

$\theta = \arccos{\frac{\vec{e} \cdot \vec{r}}{e r}}$

To place $\theta$ in the correct quadrant, it is determined first whether the spacecraft is flying away from perigee or towards perigee. If flying away from perigee ($v_r \ge 0$), then θ lies in the interval [0°, 180°), whereas if flying towards perigee ($v_r < 0$), then θ takes values on the interval [180°, 360°). In summary:

>$\begin{align}
>\theta&=\arccos{\frac{\vec{e} \cdot \vec{r}}{e r}}, (v_r \ge 0),\\
>\theta&= 360° - \arccos{\frac{\vec{e} \cdot \vec{r}}{e r}}, (v_r < 0)
>\end{align}$

"""

# ╔═╡ a5c8b870-fa82-4e5c-9934-625b755e95c8
# This function converts the state vector of a spacecraft from ECI position and velocity to COE.

function RV2COE(s, μ)
	
		R = s[1:3]
		V = s[4:6]
		
		r = norm(R) # Magnitude of position vector in geocentric equatorial plane
		v = norm(V) # Magnitude of velocity vector in geocentric equatorial plane
	
		vᵣ = dot(R, V)/r # Radial velocity component [km/s]
	
		H = cross(R, V) # Angular momentum vector [km²/s]
		h = norm(H) # magnitude of H

	a = (μ/2) * ((μ/r) - 0.5*v^2)^-1 # semi-major axis [km]

		E = cross(V, H)/μ - R/r # eccentricity vector
	e = norm(E) # eccentricity (magnitude of E)

	i = acos(H[3]/h) # inclination of orbit [rad]

		N = cross([0, 0, 1], H) # Line of nodes vector [km²/s]
		n = norm(N) # Magnitude of N

	# Right ascention of the ascending node [rad]
	if n != 0
		Ω = acos(N[1]/n)
		if N[2] < 0
			Ω = 2*pi - Ω
		end
	else
		Ω = 0
	end

	# Argument of perigee [rad]
	if n != 0
		if e > eps()
			ω = acos(dot(N, E)/(n*e))
			if E[3] < 0
				ω = 2*pi - ω
			end
		else
			ω = 0
		end
	else
		ω = 0
	end

	# True anomaly [rad]
	if e > eps()
		θ = acos(dot(E, R)/(e*r))
		if vᵣ < 0
			θ = 2*pi - θ
		end
	else
		cp = cross(N, R)
		if cp[3] >= 0
			θ = acos(dot(N, R)/(n*r))
		else
			θ = 2*pi - acos(dot(N, R)/(n*r))
		end
	end

	COE = [a, e, i, Ω, ω, θ]
end

# ╔═╡ 07c26752-1429-4cd7-bc71-f4fb3fde7b1e
# Report the COE state for the initial position given in 2.1.3:
coestate = (RV2COE(sᵢ, μₑ));

# ╔═╡ 712f4850-77cc-4e5e-b74f-bfa075e337fe
md"""
The COE state for the initial position given by $sᵢ$ is: 

$[a, e, i, \Omega, \omega, \theta] =$

[$(coestate[1]), $(coestate[2]), $(coestate[3]), $(coestate[4]), $(coestate[5]), $(coestate[6])]
	
"""

# ╔═╡ 9b7096ea-1e09-4bf3-9e0f-a34561f7ccbf
md"""
### 2.2.2
"""

# ╔═╡ 84814e56-cc4a-4f20-ab92-72100d5bf0d0
md"""
A function is written to compute the radial, transverse, normal (RTN) to Earth-centered inertial rotation matrix. If applied to a position vector in the RTN frame, the output matrix will transform that vector into the equivalent position vector in the ECI frame.

This matrix is formed such that each column contains the components of the RTN frame $[\hat{e_r}, \hat{e_t}, \hat{e_n}]$ in the ECI frame [x, y, z].

Noting that:

$\frac{\vec{r}}{r} = ê_r$

The first column is:

$\left[ \frac{\vec{r}}{r}[1], \frac{\vec{r}}{r}[2], \frac{\vec{r}}{r}[3] \right]$

Now noting that the cross product of $\vec{r}$ and $\vec{v}$ will return the out-of-plane (normal) direction:

$\frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|} = \hat{e_n}$

The third column is:

$\left[ \frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|}[1], \frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|}[2], \frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|}[3] \right]$

Finally, noting that the cross product of the first two RTN components forms the third:

$\frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|} \times \frac{\vec{r}}{r} = \hat{e_t}$

The second column is:

$\left[ (\frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|} \times \frac{\vec{r}}{r})[1], (\frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|} \times \frac{\vec{r}}{r})[2], (\frac{\vec{r} \times \vec{v}}{\|\vec{r} \times \vec{v}\|} \times \frac{\vec{r}}{r})[3] \right]$

"""

# ╔═╡ 2bc32094-7427-4e02-a76d-fc77722fc280
function RTN2ECI(s)
    r = s[1:3] # s is the Inertial state (position and velocity)
    v = s[4:6]
	
    n = cross(r, v)

    R = normalize(r)
    N = normalize(n)
    T = cross(N, R)

    R_rtn2eci = hcat(R, T, N) # Rotation matrix transforming from the RTN frame to the ECI frame

    return R_rtn2eci
end

# ╔═╡ 76cf69bb-9af7-431c-98fe-e427d40f6552
# Perform a check on the R_rtn2eci rotation matrix: inputting a starting state such that the RTN coords are already in line with the ECI coords should output the Identity matrix

RTN2ECI([1, 0, 0, 0, 1, 0])

# ╔═╡ afb89567-1d3f-4a39-bfa8-d99cffeac6dd
md"""
### 2.2.3

A ΔV impulse is applied along each direction of the RTN frame as framed at the initial condition $sᵢ$ in three seperate cases. The size of the impulse is 0.05 [km/s] for each.

The pre- and post-ΔV COE states are calculated and the difference is taken between the two to observe the affect of the ΔV on the COE of the orbit. This is done by first defining the ΔV impulse in the RTN frame. This impulse is then converted to the ECI frame via the RTN2ECI rotation matrix defined above and added to the state vector. Once added to the state vector, the state vector is converted back to COE and compared to the original state (defined in COE).
"""

# ╔═╡ 50f84de9-9092-43f9-b0c2-abe76e40342f
begin
# CASE (i): Apply a 0.05 [km/s] ΔV impulse along the radial direction
	ΔVi = [0.05, 0, 0] # ΔV impulse in RTN frame
	si = sᵢ + vcat(zeros(3), (RTN2ECI(sᵢ)*ΔVi))
	pre = RV2COE(sᵢ, μₑ) # Pre-ΔV COE state
	coei = RV2COE(si, μₑ) # Post-ΔV COE state
	Δcoei = coei - pre # Difference between pre- and post-ΔV COE states
	

# CASE (ii): Apply a 0.05 [km/s] ΔV impulse along the transverse direction
	ΔVii = [0, 0.05, 0] # ΔV impulse in RTN frame
	sii = sᵢ + vcat(zeros(3), (RTN2ECI(sᵢ)*ΔVii))
	coeii = RV2COE(sii, μₑ) # Post-ΔV COE state
	Δcoeii = coeii - pre # Difference between pre- and post-ΔV COE states

# CASE (iii): Apply a 0.05 [km/s] ΔV impulse along the normal direction
	ΔViii = [0, 0, 0.05] # ΔV impulse in RTN frame
	siii = sᵢ + vcat(zeros(3), (RTN2ECI(sᵢ)*ΔViii))
	coeiii = RV2COE(siii, μₑ) # Post-ΔV COE state
	Δcoeiii = coeiii - pre # Difference between pre- and post-ΔV COE state

[Δcoei, Δcoeii, Δcoeiii]
end;

# ╔═╡ c285ca6d-8d40-41cd-b643-be8704798336
md"
The change in COE ([$\Delta a, \Delta e, \Delta i, \Delta \Omega, \Delta \omega, \Delta \theta$]) for the radial impulse is: $(string(Δcoei))

The change in COE ([$\Delta a, \Delta e, \Delta i, \Delta \Omega, \Delta \omega, \Delta \theta$]) for the transverse impulse is: $(string(Δcoeii))

The change in COE ($\left[ \Delta a, \Delta e, \Delta i, \Delta \Omega, \Delta \omega, \Delta \theta \right]$) for the normal impulse is: $(string(Δcoeiii))

Note that the machine epsilon of the floating point type is of the order $10^{-16}$, and so the inclination does not experience any change for the radial and transverse impulses, but does for the normal impulse. This is because if the magnitudes of $v_t$ and $v_r$ remain unchanged in the impulse process, then a rigid body rotation of the orbit is produced. That is, except for its new orientation in space, the orbit remains unchanged. If the magnitudes of $v_r$ and $v_t$ change in the process, then the orbit acquires a new size and shape.

"

# ╔═╡ 61e8e74c-ccb6-46dc-a5da-46f188784f70
md"""
### 2.2.4

A grid 'search' over the true anomaly is created, whereby for each (increasing) value of θ (starting from $sᵢ$), the COE state is converted to the ECI 'RV' state, allowing the addition of the ΔV impulse using the same method as above. Once this has been applied, the RV state is converted back to COE state representation and the difference between the initial COE state for that value of true anomaly and the post-ΔV state for that value of true anomaly is computed and plotted.
"""

# ╔═╡ b2cba02d-08da-445b-ad34-8971bfcdd882
# This function from section 1.3 computes the state vector (r, v) from the classical orbital elements (COE)

function COE2RV(œ, μ)

	a, e, i, Ω, ω, θ = œ # semi-major axis, eccentricity, inclination, right ascension of the ascending node (rad), argument of perigee (rad), true anomaly (rad)

# Compute the angular momentum (h)
	
	h = sqrt(μ*a*(1-e^2)) # Note that p = a(1-e²) = h²/μ

# Compute the postion and velocity vectors in the perifocal frame
	
	rₚ = (h^2/μ) * (1/(1+e*cos(θ))) * (cos(θ)*[1, 0, 0] + sin(θ)*[0,1,0])
	vₚ = (μ/h) * (-sin(θ)*[1,0,0] + ((e+cos(θ))*[0,1,0]))

# Define the transformation matrix as the product of 3 rotation matrices about a single axis
	
	R₃Ω = [cos(Ω) sin(Ω) 0; -sin(Ω) cos(Ω) 0; 0 0 1]
	R₁i = [1 0 0; 0 cos(i) sin(i); 0 -sin(i) cos(i)]
	R₃ω = [cos(ω) sin(ω) 0; -sin(ω) cos(ω) 0; 0 0 1]

	A = transpose(R₃ω*R₁i*R₃Ω)

# Compute the state vector (r, v) in the geocentric equatorial plane

	r = A*rₚ
	v = A*vₚ

	return [r; v]
end;

# ╔═╡ 147aec25-366c-4450-a0b4-088f1fe4ead9
# Perform a check on the results by converting the COE state back to the state vector, which should yield sᵢ
COE2RV(coestate, μₑ);

# ╔═╡ 7dd51cc1-d3aa-4d94-b1d8-01205eb1aa23
md"""
The RV2COE function can be easily validated by converting back to the original state vector, which should yield sᵢ once again. It is found converting the COE state back to RV state gives:

[$(COE2RV(coestate, μₑ)[1]), $(COE2RV(coestate, μₑ)[2]), $(COE2RV(coestate, μₑ)[3]), $(COE2RV(coestate, μₑ)[4]), $(COE2RV(coestate, μₑ)[5]), $(COE2RV(coestate, μₑ)[6])]

Which is indeed equal to sᵢ (differing in rounding only).
"""

# ╔═╡ 86bf6b83-d989-47e7-a758-64b458d826a8
begin
	θ = LinRange(0.01, 2*pi, 1000)
	grid = zeros(1000, 18)

	for i = 1:1000
		PRE = vcat(RV2COE(sᵢ, μₑ)[1:5], θ[i])
		
		COEi = RV2COE((COE2RV(PRE, μₑ) + vcat(zeros(3), (RTN2ECI(COE2RV(PRE, μₑ))*ΔVi))), μₑ)
		COEii = RV2COE((COE2RV(PRE, μₑ) + vcat(zeros(3), (RTN2ECI(COE2RV(PRE, μₑ))*ΔVii))), μₑ)
		COEiii = RV2COE((COE2RV(PRE, μₑ) + vcat(zeros(3), (RTN2ECI(COE2RV(PRE, μₑ))*ΔViii))), μₑ)
		
		ΔCOEi = COEi - PRE
		ΔCOEii = COEii - PRE
		ΔCOEiii = COEiii - PRE
		
		grid[i, :] .= vcat(ΔCOEi, ΔCOEii, ΔCOEiii)
	end
end

# ╔═╡ 18254000-21c0-407d-8d31-080321b51c6d
begin
# Pull out solutions for plotting:
	ar = grid[:, 1]
	at = grid[:, 7]
	an = grid[:, 13]
	
	er = grid[:, 2]
	et = grid[:, 8]
	en = grid[:, 14]
	
	ir = grid[:, 3]
	it = grid[:, 9]
	i_n = grid[:, 15]
	
	Ωr = grid[:, 4]
	Ωt = grid[:, 10]
	Ωn = grid[:, 16]
	
	ωr = grid[:, 5]
	ωt = grid[:, 11]
	ωn = grid[:, 17]
	
	θr = grid[:, 6]
	θt = grid[:, 12]
	θn = grid[:, 18]
end;

# ╔═╡ 3915c82f-f36b-4ca7-9357-c799f25bc0a1
md"""
Noting that perigee of the orbit lies at $\theta = 0$ (or $\theta = 2 \pi$), and the apogee lies at $\theta = \pi$, the maximal change in semi-major axis occurs at perigee when the transverse impulse is applied. The minimum change occurs at apogee. This agrees with the fact that the energy change of the orbit is maximised:

1) When the orbit velocity is maximal (ie. at perigee)
2) When the impulse is applied tangentially

These two points are clear upon inspection of the change in orbital energy formula:

$\Delta E = \frac{1}{2} \Delta v^2 + v_0 \Delta v \cos{\alpha}$

Where $\alpha$ is the angle between the orignial velocity vector $v_0$ and $\Delta v$.

The purely normal impulse has no affect on the size of the orbit as it affects solely the orientation of the orbit as it points out-of-plane.
"""

# ╔═╡ 22643b1c-8b92-4e46-b071-7f4c99019949
begin
	θd = rad2deg.(θ)
	
	plot_a = plot(
		θd,
		ar,
		label = "radial",
		title = "Change in semi-major axis due to delta-V impulse for range of true anomaly",
		xlabel = "True Anomaly [deg]",
		ylabel = "Change in semi-major axis [km]",
		titlefontsize = 9,
		legend = :right,
		legendcolumns = 1,
	)
	plot!(
		θd,
		at,
		label = "transverse"
	)
	plot!(
		θd,
		an,
		label = "normal"
	)
end

# ╔═╡ 2b296a5c-6a47-459c-8b56-8100850c59ca
md"""
Due to the same reason as discussed for the semi-major axis, the maximal change in eccentricity also occurs with a transverse burn at perigee, and no change is incurred with a normal impulse.
"""

# ╔═╡ 903a72db-3357-4018-ba28-2dfb5905714b
begin
	plot_e = plot(
		θd,
		er,
		label = "radial",
		title = "Change in eccentricity due to delta-V impulse for range of true anomaly",
		xlabel = "True Anomaly [deg]",
		ylabel = "Change in eccentricity",
		titlefontsize = 9,
		legend = :bottomright,
		legendcolumns = 1,
	)
	plot!(
		θd,
		et,
		label = "transverse"
	)
	plot!(
		θd,
		en,
		label = "normal"
	)
end

# ╔═╡ fd606282-1ff2-4d78-a120-2ac660b219fa
begin
	plot_i = plot(
		θd,
		rad2deg.(ir),
		label = "radial",
		title = "Change in inclination due to delta-V impulse for range of true anomaly",
		xlabel = "True Anomaly [deg]",
		ylabel = "Change in inclination [deg]",
		titlefontsize = 9,
		legend = :bottomright,
		legendcolumns = 1,
	)
	plot!(
		θd,
		rad2deg.(it),
		label = "transverse"
	)
	plot!(
		θd,
		rad2deg.(i_n),
		label = "normal"
	)
end

# ╔═╡ 0345834c-1818-4023-9de8-04f732269dba
begin
	plot_Ω = plot(
		θd,
		rad2deg.(Ωr),
		label = "radial",
		title = "Change in Ω due to delta-V impulse for range of true anomaly",
		xlabel = "True Anomaly [deg]",
		ylabel = "Change in Right Ascension of Ascending Node [deg]",
		titlefontsize = 10,
		ylabelfontsize = 8,
		xlabelfontsize = 10,
		legend = :bottomright,
		legendcolumns = 1,
	)
	plot!(
		θd,
		rad2deg.(Ωt),
		label = "transverse"
	)
	plot!(
		θd,
		rad2deg.(Ωn),
		label = "normal"
	)
end

# ╔═╡ 07692bc0-e289-4ac3-af4e-2ee7ec29b6a2
begin
	plot_ω = plot(
		θd,
		rad2deg.(ωr),
		label = "radial",
		title = "Change in ω due to delta-V impulse for range of true anomaly",
		xlabel = "True Anomaly [deg]",
		ylabel = "Change in ω [deg]",
		titlefontsize = 10,
		ylabelfontsize = 10,
		xlabelfontsize = 10,
		legend = :topright,
		legendcolumns = 1,
	)
	plot!(
		θd,
		rad2deg.(ωt),
		label = "transverse"
	)
	plot!(
		θd,
		rad2deg.(ωn),
		label = "normal"
	)
end

# ╔═╡ 861c4fb6-13cf-43db-b8f2-4cff0142dcd2
begin
	plot_θ = plot(
		θd,
		rad2deg.(θr),
		label = "radial",
		title = "Change in true anomaly due to delta-V impulse for range of true anomaly",
		xlabel = "True Anomaly [deg]",
		ylabel = "Change in θ [deg]",
		titlefontsize = 10,
		ylabelfontsize = 10,
		xlabelfontsize = 10,
		legend = :topright,
		legendcolumns = 1,
	)
	plot!(
		θd,
		rad2deg.(θt),
		label = "transverse"
	)
	plot!(
		θd,
		rad2deg.(θn),
		label = "normal"
	)
end;

# ╔═╡ fa897a98-3304-417f-9074-7c84a1c20356
md"""
### 2.3.1
"""

# ╔═╡ de757199-f956-4e69-bdc0-7910b971ef34
md"""
The average rate of precession of the node line, and hence, the orbital plane is given by:

$\dot{\Omega} = - \left[ \frac{3}{2} \frac{\sqrt{\mu} J_2 R^2}{(1-e^2) a^{\frac{7}{2}}} \right] \cos{i}$

Where $R$ and $μ$ are the radius and gravitational parameter of the Earth, $a$ and $e$ are the semimajor axis and eccentricity of the orbit, and $i$ is the orbit’s inclination. Note that if $0 ≤ i < 90°$, then $\dot{\Omega} < 0$, i.e. for prograde orbits, the node line drifts westward. Conversely, If $90\degree < i ≤ 180\degree$, we see that $\dot{\Omega} > 0 \implies$ the node line of retrograde orbits advance eastward.

For a sun-synchronous orbit, the orbital plane must rotate in inertial space with the angular velocity of the earth in its orbit around the sun, which is 360° per 365.26 days $= 0.9856\degree$ per day. Observe that this imposes the requirement that the orbital plane precesses eastward i.e. $\dot{\Omega} > 0$, and therefore, the orbit must be retrograde with $i > 90\degree$. The inclination at which the orbit is sun-synchronous may then be obtained upon rearranging the formula for $\dot{\Omega}$ for $i$.
"""

# ╔═╡ 212b15dd-5ba1-49f8-8fbc-c28e8dcfa69a
begin
		R = sᵢ[1:3]
		V = sᵢ[4:6]
		r = norm(R) # Magnitude of position vector in geocentric equatorial plane
		v = norm(V) # Magnitude of velocity vector in geocentric equatorial plane
		H = cross(R, V) # Angular momentum vector [km²/s]
		h = norm(H) # magnitude of H

	a = (μₑ/2) * ((μₑ/r) - 0.5*v^2)^-1 # semi-major axis [km]
		E = cross(V, H)/μₑ - R/r # eccentricity vector
	e = norm(E) # eccentricity (magnitude of E)

	dΩ = 1.991e-7 # Ascending node advance rate [rad/s]

	i = acos(dΩ/(- (3/2 * (sqrt(μₑ) * J₂ * Rₑ^2)/((1-e^2) * a^(7/2))))) # Inclination required for sun-synchronous orbit

	i_deg = rad2deg(i)
end;

# ╔═╡ 82cee6e0-36f5-4ff4-9517-3c4ce95c4fa5
md"""
The inclination required to place the satellite in sun-synchronous orbit is: $(i_deg) [deg].

To check this result satisfies the requirements of sun-synchronous orbits: $\cos(i) > 0 \implies i > 90\degree$ which is indeed satisfied by the value calculated.
"""

# ╔═╡ 43ff924f-8b74-4947-ae4d-5b6c28b58eba
md"""
### 2.3.2
"""

# ╔═╡ 8b5a4da4-85c7-41f7-8568-cd89d6ab2eea
md"""
The velocity of the spacecraft may be determined at the ascending and descending nodes using the equation of the orbit and the Vis-Viva equation. First the location of the ascending node is determined from the argument of perigee as follows:

$\theta_{ascend} = 2 \pi - \omega$
$\theta_{descend} = \pi - \omega$

The equation of the orbit may then be applied at each node as:

$r = \frac{\frac{h^2}{\mu}}{1+e \cos{\theta}}$

This gives the magnitude of the radial position of each node, which may then be fed into the Vis-Viva equation:

$\frac{1}{2} v^2 - \frac{\mu}{r} = - \frac{\mu}{2 a}$

and solved for $v$:

$v = \left[ 2 \left( \frac{\mu}{r} - \frac{\mu}{2 a} \right) \right] ^{\frac{1}{2}}$


"""

# ╔═╡ 0ba4cf06-cb37-44a9-8963-92fa3c660dc4
begin
	orbit = RV2COE(sᵢ, μₑ)
	omega = orbit[5]
	theta_ascend = 2*pi-omega
	theta_descend = pi-omega
	
r_asc = (h^2/μₑ)/(1+e*cos(theta_ascend))
r_des = (h^2/μₑ)/(1+e*cos(theta_descend))

v_asc = sqrt(2*(μₑ/r_asc - μₑ/(2*a)))
v_des = sqrt(2*(μₑ/r_des - μₑ/(2*a)))
end;

# ╔═╡ f3216d7f-6978-4123-8e72-f28b859dea26
md"""
The value of True Anomaly at the ascending and descending node are $(rad2deg(theta_ascend)) [deg] and $(rad2deg(theta_descend)) [deg] respectively. 

The maximum inclination change occurs at the ascending and desceding nodes. Note these points are those at which the argument of latitude $u = \theta + \omega$ is 0 or 180°.

The reason for this is that the line of nodes is what the orbit tilt, or inclination, rotates about. Applying the delta-V impulse at these points ensures that no delta-V component is spent rotating the orbit around the vertical polar axis, and a pure change in the inclination is incurred, thus it is most efficient here. Applying the delta-V at any other point on the orbit would mean some of the out-of-plane impulse would create a velocity component that rotates the line of nodes with respect to the ECI frame, which changes the Ω of the orbit aswell as the inclination. Indeed, if the delta-V in the $\hat{e_n}$ direction were to be applied at the location in-line with the north or south poles of the Earth, then no inclination change would be incurred, and the orbit rotation would occur purely around the vertical polar axis, changing the Ω of the orbit only. In summary, the maximum change of Ω is incurred at the maximum and minimum latitudes, while the minimum (zero) effect is observed on the node line at the ascending and descending nodes. Conversely, the maximum change of $i$ is incurred at the ascending and descending nodes, while the minimum (zero) effect is observed at the maximum and minimum latitudes.
"""

# ╔═╡ 2a5ea354-efa2-45e5-9f60-9c4e3d77540a
begin # Check against a different method
	b1 = COE2RV(vcat(orbit[1:5], theta_ascend), μₑ)
	b2 = COE2RV(vcat(orbit[1:5], theta_descend), μₑ)
	v1 = norm(b1[4:6])
	v2 = norm(b2[4:6])

	# Return 'true' if the difference of the velocities is less than floating pt error
	[v1, v2] - [v_asc, v_des] <= [1e-16, 1e-16]
end;

# ╔═╡ 597c6890-f1d3-4f2d-9031-9c834b7d9968
# perform a check on the selection of the ascending and descending nodes
begin
	coe_ascend = replace(orbit, orbit[6] => theta_ascend)
	coe_descend = replace(orbit, orbit[6] => theta_descend)

# z component of v should be positive at the ascending node and negative at the descending node in the ECI frame
	
	rva = COE2RV(coe_ascend, μₑ)
	rvd = COE2RV(coe_descend, μₑ)
	[rva[6], rvd[6]]
end;

# ╔═╡ 2271c992-8f4a-4ad8-a04a-111404ed2576
md"""
The value of True Anomaly at the ascending and descending node are $(rad2deg(theta_ascend)) [deg] and $(rad2deg(theta_descend)) [deg] respectively.

Values of the velocity of the spacecraft at the ascending and descending node are $(v_asc) [km/s] and $(v_des) [km/s] respectively.

A quick check on the selection of the ascending and descending nodes is to verify that the z component of the velocity vector in the ECI frame is positive at the ascending node and negative at the descending node. It is found that at the ascending node, the z-component is: $(rva[6]) [km/s] and at the descending node, it is: $(rvd[6]) [km/s]. This passes the check.
"""

# ╔═╡ 6a21606e-8d0d-4502-97c4-9895b39f9d21
md"""
### 2.3.3
The ΔV required to change the inclination of the spacecrafts orbit from initial inclination $(orbit[3]) to the sun-synchronous orbit inclination $(i), is given by the following formula for delta-v for a pure rotation of the velocity vector through the inclination angle. 

$\Delta V = 2 v_θ \sin{\frac{\Delta i}{2}}$

As such, the best of the two above nodes to execute an inclination change maneuver in terms of ΔV efficiency, will be the one at which the transverse component of the velocity $v_{\theta}$ is lowest. This corresponds to the node at which the magnitude of $\vec{v}$ is lowest.

The values of the velocity of the spacecraft at the ascending and descending node were $(v_asc) [km/s] and $(v_des) [km/s] respectively, and so the best node to perform the inclination change manoeuvre is the descending node.

In order to compute $v_{\theta}$, the velocity vector at the descending node must be converted to RTN coordinates via multiplication with the transpose of the previously computed RTN2ECI matrix. The second component of the velocity in RTN coordinates is the desired $v_{\theta}$. This may then be used in $\Delta V = 2 v_θ \sin{\frac{\Delta i}{2}}$ to compute the required ΔV.
"""

# ╔═╡ 3f613597-abc4-4b38-ae03-238856fc89a3
begin
	Δi = abs(orbit[3]-i);

	rtn_des = (RTN2ECI(rvd))' * rvd[4:6];
	vₜ = rtn_des[2];
	
ΔV = 2*vₜ*sin(Δi/2);
percent = ΔV/(norm(rvd[4:6]))*100;
	
end;

# ╔═╡ 499dcecc-8c23-4570-83b5-b776517f665a
md"""
The ΔV required at the descending node is $(ΔV) [km/s]. As a percentage of the original velocity at the descending node, this is $(percent)%.
"""

# ╔═╡ 1f5b5feb-6ccd-4e6c-8979-1c688409924d
md"""
### 2.3.4
Given that the Sun-synchronous orbit computed above was computed by a pure rotation about the line of nodes, all the orbital elements remain the same except for the inclination. As such, the $a$ and $e$ of the original orbit remain the same for the sun-synchronized one, and thus these values may be used to calculate $rₐ$ , the apogee distance.

$rₐ = a(1+e)$

A new set of orbits may then be defined with the same rₐ as the original sun-synchronous orbit, but a different set of values for rₚ (and therefore different semi-major axes and eccentricities).

The Vis-Viva equation may then be used to calculate the new value of $v$ at apogee for each orbit in the set.

The ∆V required to lower the perigee to a range of altitudes from 60 km to 210
km above the Earth’s surface is then computed from the difference between the 'new' values of $v$ at apogee and the original value of $v$ at apogee.

"""

# ╔═╡ 7429fc5e-7275-4ae1-876b-0fa83a5714c0
begin
# Note a and e of sun-synchronous orbit are unchanged from the original orbit.

	rₐ = a*(1+e) # Radius at apogee for sun-synchronous orbit

	alt = LinRange(60, 210, 16)
	rₚ = Rₑ .+ alt

# Define a new set of orbits with the same rₐ as the original sun-synchronous orbit, but a different set of values for rₚ (and therefore different semi-major axis and eccentricity).

	â = (rₚ .+ rₐ)/2

# Use the Vis-Viva equation to calculate the new value of v at apogee for each orbit in the set.

	v̂ = sqrt.(2 .* ((μₑ/rₐ) .- (μₑ./(2 .* â))))

# Compute the initial value of v at apogee

	vᵢ = sqrt(2 * ((μₑ/rₐ) - (μₑ/(2 * a))))

# Compute the ΔV required to change from the original orbit to each new orbit

	ΔV̂ = abs.(v̂ .- vᵢ)

# Plot the ∆V required to lower the perigee to a range of 16 altitudes from 60 km to 210 km above the Earth’s surface

	plot(
			alt,
			ΔV̂,
			legend = false,
			titlefontsize = 12,
			title = "ΔV Required for Perigee Reduction to given Altitude",
			xlabel = "Altitude [km]",
			xlabelfontsize = 10,
			ylabel = "ΔV [km/s]",
			ylabelfontsize = 10
		)
	
end

# ╔═╡ 04f92c6d-b580-44ea-a9c8-41e13fe02d25
md"""
Since the velocities at apogee of the 'new' orbits are less than those of the original orbit, the ΔV impulses performed to transfer to the new orbits are accomplished by retrofires. That is, the thrust of the maneuvering spacecraft is directed opposite to the flight direction to act as a brake on the motion. Since ΔV represents the same propellent expenditure regardless of the direction the thruster is aimed, the magnitude of the ΔV is plotted.

The perigee of the original orbit was $(a*(1-e)) [km] which is equivalent to an altitude above the Earth's surface of $(a*(1-e)-Rₑ) [km]. As such, it is expected that the ΔV required to reduce the perigee altitude of the orbit to altitudes closer to the original altitude would be less, as less energy need be expended to change the orbit. This linear relationship between ΔV and perigee altitude is reflected in the results, which show a clear linear trend of decreasing ΔV with altitude increasing to the value of the original perigee altitude.


"""

# ╔═╡ 00000000-0000-0000-0000-000000000001
PLUTO_PROJECT_TOML_CONTENTS = """
[deps]
LinearAlgebra = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"
Plots = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
PlutoUI = "7f904dfe-b85e-4ff6-b463-dae2292396a8"

[compat]
Plots = "~1.40.4"
PlutoUI = "~0.7.59"
"""

# ╔═╡ 00000000-0000-0000-0000-000000000002
PLUTO_MANIFEST_TOML_CONTENTS = """
# This file is machine-generated - editing it directly is not advised

julia_version = "1.10.0"
manifest_format = "2.0"
project_hash = "634710e806d44fb13e8e9c0131ce0d1dfcd55dbd"

[[deps.AbstractPlutoDingetjes]]
deps = ["Pkg"]
git-tree-sha1 = "297b6b41b66ac7cbbebb4a740844310db9fd7b8c"
uuid = "6e696c72-6542-2067-7265-42206c756150"
version = "1.3.1"

[[deps.ArgTools]]
uuid = "0dad84c5-d112-42e6-8d28-ef12dabb789f"
version = "1.1.1"

[[deps.Artifacts]]
uuid = "56f22d72-fd6d-98f1-02f0-08ddc0907c33"

[[deps.Base64]]
uuid = "2a0f44e3-6c83-55bd-87e4-b1978d98bd5f"

[[deps.BitFlags]]
git-tree-sha1 = "2dc09997850d68179b69dafb58ae806167a32b1b"
uuid = "d1d4a3ce-64b1-5f1a-9ba4-7e7e69966f35"
version = "0.1.8"

[[deps.Bzip2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "9e2a6b69137e6969bab0152632dcb3bc108c8bdd"
uuid = "6e34b625-4abd-537c-b88f-471c36dfa7a0"
version = "1.0.8+1"

[[deps.Cairo_jll]]
deps = ["Artifacts", "Bzip2_jll", "CompilerSupportLibraries_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "JLLWrappers", "LZO_jll", "Libdl", "Pixman_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a4c43f59baa34011e303e76f5c8c91bf58415aaf"
uuid = "83423d85-b0ee-5818-9007-b63ccbeb887a"
version = "1.18.0+1"

[[deps.CodecZlib]]
deps = ["TranscodingStreams", "Zlib_jll"]
git-tree-sha1 = "59939d8a997469ee05c4b4944560a820f9ba0d73"
uuid = "944b1d66-785c-5afd-91f1-9de20f533193"
version = "0.7.4"

[[deps.ColorSchemes]]
deps = ["ColorTypes", "ColorVectorSpace", "Colors", "FixedPointNumbers", "PrecompileTools", "Random"]
git-tree-sha1 = "67c1f244b991cad9b0aa4b7540fb758c2488b129"
uuid = "35d6a980-a343-548e-a6ea-1d62b119f2f4"
version = "3.24.0"

[[deps.ColorTypes]]
deps = ["FixedPointNumbers", "Random"]
git-tree-sha1 = "b10d0b65641d57b8b4d5e234446582de5047050d"
uuid = "3da002f7-5984-5a60-b8a6-cbb66c0b333f"
version = "0.11.5"

[[deps.ColorVectorSpace]]
deps = ["ColorTypes", "FixedPointNumbers", "LinearAlgebra", "Requires", "Statistics", "TensorCore"]
git-tree-sha1 = "a1f44953f2382ebb937d60dafbe2deea4bd23249"
uuid = "c3611d14-8923-5661-9e6a-0046d554d3a4"
version = "0.10.0"

    [deps.ColorVectorSpace.extensions]
    SpecialFunctionsExt = "SpecialFunctions"

    [deps.ColorVectorSpace.weakdeps]
    SpecialFunctions = "276daf66-3868-5448-9aa4-cd146d93841b"

[[deps.Colors]]
deps = ["ColorTypes", "FixedPointNumbers", "Reexport"]
git-tree-sha1 = "fc08e5930ee9a4e03f84bfb5211cb54e7769758a"
uuid = "5ae59095-9a9b-59fe-a467-6f913c188581"
version = "0.12.10"

[[deps.Compat]]
deps = ["TOML", "UUIDs"]
git-tree-sha1 = "c955881e3c981181362ae4088b35995446298b80"
uuid = "34da2185-b29b-5c13-b0c7-acf172513d20"
version = "4.14.0"
weakdeps = ["Dates", "LinearAlgebra"]

    [deps.Compat.extensions]
    CompatLinearAlgebraExt = "LinearAlgebra"

[[deps.CompilerSupportLibraries_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "e66e0078-7015-5450-92f7-15fbd957f2ae"
version = "1.0.5+1"

[[deps.ConcurrentUtilities]]
deps = ["Serialization", "Sockets"]
git-tree-sha1 = "6cbbd4d241d7e6579ab354737f4dd95ca43946e1"
uuid = "f0e56b4a-5159-44fe-b623-3e5288b988bb"
version = "2.4.1"

[[deps.Contour]]
git-tree-sha1 = "439e35b0b36e2e5881738abc8857bd92ad6ff9a8"
uuid = "d38c429a-6771-53c6-b99e-75d170b6e991"
version = "0.6.3"

[[deps.DataAPI]]
git-tree-sha1 = "abe83f3a2f1b857aac70ef8b269080af17764bbe"
uuid = "9a962f9c-6df0-11e9-0e5d-c546b8b5ee8a"
version = "1.16.0"

[[deps.DataStructures]]
deps = ["Compat", "InteractiveUtils", "OrderedCollections"]
git-tree-sha1 = "1d0a14036acb104d9e89698bd408f63ab58cdc82"
uuid = "864edb3b-99cc-5e75-8d2d-829cb0a9cfe8"
version = "0.18.20"

[[deps.Dates]]
deps = ["Printf"]
uuid = "ade2ca70-3891-5945-98fb-dc099432e06a"

[[deps.DelimitedFiles]]
deps = ["Mmap"]
git-tree-sha1 = "9e2f36d3c96a820c678f2f1f1782582fcf685bae"
uuid = "8bb1440f-4735-579b-a4ab-409b98df4dab"
version = "1.9.1"

[[deps.DocStringExtensions]]
deps = ["LibGit2"]
git-tree-sha1 = "2fb1e02f2b635d0845df5d7c167fec4dd739b00d"
uuid = "ffbed154-4ef7-542d-bbb7-c09d3a79fcae"
version = "0.9.3"

[[deps.Downloads]]
deps = ["ArgTools", "FileWatching", "LibCURL", "NetworkOptions"]
uuid = "f43a241f-c20a-4ad4-852c-f6b1247861c6"
version = "1.6.0"

[[deps.EpollShim_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8e9441ee83492030ace98f9789a654a6d0b1f643"
uuid = "2702e6a9-849d-5ed8-8c21-79e8b8f9ee43"
version = "0.0.20230411+0"

[[deps.ExceptionUnwrapping]]
deps = ["Test"]
git-tree-sha1 = "dcb08a0d93ec0b1cdc4af184b26b591e9695423a"
uuid = "460bff9d-24e4-43bc-9d9f-a8973cb893f4"
version = "0.1.10"

[[deps.Expat_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "4558ab818dcceaab612d1bb8c19cee87eda2b83c"
uuid = "2e619515-83b5-522b-bb60-26c02a35a201"
version = "2.5.0+0"

[[deps.FFMPEG]]
deps = ["FFMPEG_jll"]
git-tree-sha1 = "b57e3acbe22f8484b4b5ff66a7499717fe1a9cc8"
uuid = "c87230d0-a227-11e9-1b43-d7ebe4e7570a"
version = "0.4.1"

[[deps.FFMPEG_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "JLLWrappers", "LAME_jll", "Libdl", "Ogg_jll", "OpenSSL_jll", "Opus_jll", "PCRE2_jll", "Zlib_jll", "libaom_jll", "libass_jll", "libfdk_aac_jll", "libvorbis_jll", "x264_jll", "x265_jll"]
git-tree-sha1 = "466d45dc38e15794ec7d5d63ec03d776a9aff36e"
uuid = "b22a6f82-2f65-5046-a5b2-351ab43fb4e5"
version = "4.4.4+1"

[[deps.FileWatching]]
uuid = "7b1f6079-737a-58dc-b8bc-7a2ca5c1b5ee"

[[deps.FixedPointNumbers]]
deps = ["Statistics"]
git-tree-sha1 = "335bfdceacc84c5cdf16aadc768aa5ddfc5383cc"
uuid = "53c48c17-4a7d-5ca2-90c5-79b7896eea93"
version = "0.8.4"

[[deps.Fontconfig_jll]]
deps = ["Artifacts", "Bzip2_jll", "Expat_jll", "FreeType2_jll", "JLLWrappers", "Libdl", "Libuuid_jll", "Pkg", "Zlib_jll"]
git-tree-sha1 = "21efd19106a55620a188615da6d3d06cd7f6ee03"
uuid = "a3f928ae-7b40-5064-980b-68af3947d34b"
version = "2.13.93+0"

[[deps.Format]]
git-tree-sha1 = "9c68794ef81b08086aeb32eeaf33531668d5f5fc"
uuid = "1fa38f19-a742-5d3f-a2b9-30dd87b9d5f8"
version = "1.3.7"

[[deps.FreeType2_jll]]
deps = ["Artifacts", "Bzip2_jll", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d8db6a5a2fe1381c1ea4ef2cab7c69c2de7f9ea0"
uuid = "d7e528f0-a631-5988-bf34-fe36492bcfd7"
version = "2.13.1+0"

[[deps.FriBidi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "aa31987c2ba8704e23c6c8ba8a4f769d5d7e4f91"
uuid = "559328eb-81f9-559d-9380-de523a88c83c"
version = "1.0.10+0"

[[deps.GLFW_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libglvnd_jll", "Xorg_libXcursor_jll", "Xorg_libXi_jll", "Xorg_libXinerama_jll", "Xorg_libXrandr_jll"]
git-tree-sha1 = "ff38ba61beff76b8f4acad8ab0c97ef73bb670cb"
uuid = "0656b61e-2033-5cc2-a64a-77c0f6c09b89"
version = "3.3.9+0"

[[deps.GR]]
deps = ["Artifacts", "Base64", "DelimitedFiles", "Downloads", "GR_jll", "HTTP", "JSON", "Libdl", "LinearAlgebra", "Pkg", "Preferences", "Printf", "Random", "Serialization", "Sockets", "TOML", "Tar", "Test", "UUIDs", "p7zip_jll"]
git-tree-sha1 = "3437ade7073682993e092ca570ad68a2aba26983"
uuid = "28b8d3ca-fb5f-59d9-8090-bfdbd6d07a71"
version = "0.73.3"

[[deps.GR_jll]]
deps = ["Artifacts", "Bzip2_jll", "Cairo_jll", "FFMPEG_jll", "Fontconfig_jll", "FreeType2_jll", "GLFW_jll", "JLLWrappers", "JpegTurbo_jll", "Libdl", "Libtiff_jll", "Pixman_jll", "Qt6Base_jll", "Zlib_jll", "libpng_jll"]
git-tree-sha1 = "a96d5c713e6aa28c242b0d25c1347e258d6541ab"
uuid = "d2c73de3-f751-5644-a686-071e5b155ba9"
version = "0.73.3+0"

[[deps.Gettext_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "Libdl", "Libiconv_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "9b02998aba7bf074d14de89f9d37ca24a1a0b046"
uuid = "78b55507-aeef-58d4-861c-77aaff3498b1"
version = "0.21.0+0"

[[deps.Glib_jll]]
deps = ["Artifacts", "Gettext_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Libiconv_jll", "Libmount_jll", "PCRE2_jll", "Zlib_jll"]
git-tree-sha1 = "359a1ba2e320790ddbe4ee8b4d54a305c0ea2aff"
uuid = "7746bdde-850d-59dc-9ae8-88ece973131d"
version = "2.80.0+0"

[[deps.Graphite2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "344bf40dcab1073aca04aa0df4fb092f920e4011"
uuid = "3b182d85-2403-5c21-9c21-1e1f0cc25472"
version = "1.3.14+0"

[[deps.Grisu]]
git-tree-sha1 = "53bb909d1151e57e2484c3d1b53e19552b887fb2"
uuid = "42e2da0e-8278-4e71-bc24-59509adca0fe"
version = "1.0.2"

[[deps.HTTP]]
deps = ["Base64", "CodecZlib", "ConcurrentUtilities", "Dates", "ExceptionUnwrapping", "Logging", "LoggingExtras", "MbedTLS", "NetworkOptions", "OpenSSL", "Random", "SimpleBufferStream", "Sockets", "URIs", "UUIDs"]
git-tree-sha1 = "2c3ec1f90bb4a8f7beafb0cffea8a4c3f4e636ab"
uuid = "cd3eb016-35fb-5094-929b-558a96fad6f3"
version = "1.10.6"

[[deps.HarfBuzz_jll]]
deps = ["Artifacts", "Cairo_jll", "Fontconfig_jll", "FreeType2_jll", "Glib_jll", "Graphite2_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg"]
git-tree-sha1 = "129acf094d168394e80ee1dc4bc06ec835e510a3"
uuid = "2e76f6c2-a576-52d4-95c1-20adfe4de566"
version = "2.8.1+1"

[[deps.Hyperscript]]
deps = ["Test"]
git-tree-sha1 = "179267cfa5e712760cd43dcae385d7ea90cc25a4"
uuid = "47d2ed2b-36de-50cf-bf87-49c2cf4b8b91"
version = "0.0.5"

[[deps.HypertextLiteral]]
deps = ["Tricks"]
git-tree-sha1 = "7134810b1afce04bbc1045ca1985fbe81ce17653"
uuid = "ac1192a8-f4b3-4bfe-ba22-af5b92cd3ab2"
version = "0.9.5"

[[deps.IOCapture]]
deps = ["Logging", "Random"]
git-tree-sha1 = "8b72179abc660bfab5e28472e019392b97d0985c"
uuid = "b5f81e59-6552-4d32-b1f0-c071b021bf89"
version = "0.2.4"

[[deps.InteractiveUtils]]
deps = ["Markdown"]
uuid = "b77e0a4c-d291-57a0-90e8-8db25a27a240"

[[deps.IrrationalConstants]]
git-tree-sha1 = "630b497eafcc20001bba38a4651b327dcfc491d2"
uuid = "92d709cd-6900-40b7-9082-c6be49f344b6"
version = "0.2.2"

[[deps.JLFzf]]
deps = ["Pipe", "REPL", "Random", "fzf_jll"]
git-tree-sha1 = "a53ebe394b71470c7f97c2e7e170d51df21b17af"
uuid = "1019f520-868f-41f5-a6de-eb00f4b6a39c"
version = "0.1.7"

[[deps.JLLWrappers]]
deps = ["Artifacts", "Preferences"]
git-tree-sha1 = "7e5d6779a1e09a36db2a7b6cff50942a0a7d0fca"
uuid = "692b3bcd-3c85-4b1f-b108-f13ce0eb3210"
version = "1.5.0"

[[deps.JSON]]
deps = ["Dates", "Mmap", "Parsers", "Unicode"]
git-tree-sha1 = "31e996f0a15c7b280ba9f76636b3ff9e2ae58c9a"
uuid = "682c06a0-de6a-54ab-a142-c8b1cf79cde6"
version = "0.21.4"

[[deps.JpegTurbo_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3336abae9a713d2210bb57ab484b1e065edd7d23"
uuid = "aacddb02-875f-59d6-b918-886e6ef4fbf8"
version = "3.0.2+0"

[[deps.LAME_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "f6250b16881adf048549549fba48b1161acdac8c"
uuid = "c1c5ebd0-6772-5130-a774-d5fcae4a789d"
version = "3.100.1+0"

[[deps.LERC_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "bf36f528eec6634efc60d7ec062008f171071434"
uuid = "88015f11-f218-50d7-93a8-a6af411a945d"
version = "3.0.0+1"

[[deps.LLVMOpenMP_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "d986ce2d884d49126836ea94ed5bfb0f12679713"
uuid = "1d63c593-3942-5779-bab2-d838dc0a180e"
version = "15.0.7+0"

[[deps.LZO_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "e5b909bcf985c5e2605737d2ce278ed791b89be6"
uuid = "dd4b983a-f0e5-5f8d-a1b7-129d4a5fb1ac"
version = "2.10.1+0"

[[deps.LaTeXStrings]]
git-tree-sha1 = "50901ebc375ed41dbf8058da26f9de442febbbec"
uuid = "b964fa9f-0449-5b57-a5c2-d3ea65f4040f"
version = "1.3.1"

[[deps.Latexify]]
deps = ["Format", "InteractiveUtils", "LaTeXStrings", "MacroTools", "Markdown", "OrderedCollections", "Requires"]
git-tree-sha1 = "e0b5cd21dc1b44ec6e64f351976f961e6f31d6c4"
uuid = "23fbe1c1-3f47-55db-b15f-69d7ec21a316"
version = "0.16.3"

    [deps.Latexify.extensions]
    DataFramesExt = "DataFrames"
    SymEngineExt = "SymEngine"

    [deps.Latexify.weakdeps]
    DataFrames = "a93c6f00-e57d-5684-b7b6-d8193f3e46c0"
    SymEngine = "123dc426-2d89-5057-bbad-38513e3affd8"

[[deps.LibCURL]]
deps = ["LibCURL_jll", "MozillaCACerts_jll"]
uuid = "b27032c2-a3e7-50c8-80cd-2d36dbcbfd21"
version = "0.6.4"

[[deps.LibCURL_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll", "Zlib_jll", "nghttp2_jll"]
uuid = "deac9b47-8bc7-5906-a0fe-35ac56dc84c0"
version = "8.4.0+0"

[[deps.LibGit2]]
deps = ["Base64", "LibGit2_jll", "NetworkOptions", "Printf", "SHA"]
uuid = "76f85450-5226-5b5a-8eaa-529ad045b433"

[[deps.LibGit2_jll]]
deps = ["Artifacts", "LibSSH2_jll", "Libdl", "MbedTLS_jll"]
uuid = "e37daf67-58a4-590a-8e99-b0245dd2ffc5"
version = "1.6.4+0"

[[deps.LibSSH2_jll]]
deps = ["Artifacts", "Libdl", "MbedTLS_jll"]
uuid = "29816b5a-b9ab-546f-933c-edad1886dfa8"
version = "1.11.0+1"

[[deps.Libdl]]
uuid = "8f399da3-3557-5675-b5ff-fb832c97cbdb"

[[deps.Libffi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "0b4a5d71f3e5200a7dff793393e09dfc2d874290"
uuid = "e9f186c6-92d2-5b65-8a66-fee21dc1b490"
version = "3.2.2+1"

[[deps.Libgcrypt_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgpg_error_jll", "Pkg"]
git-tree-sha1 = "64613c82a59c120435c067c2b809fc61cf5166ae"
uuid = "d4300ac3-e22c-5743-9152-c294e39db1e4"
version = "1.8.7+0"

[[deps.Libglvnd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll", "Xorg_libXext_jll"]
git-tree-sha1 = "6f73d1dd803986947b2c750138528a999a6c7733"
uuid = "7e76a0d4-f3c7-5321-8279-8d96eeed0f29"
version = "1.6.0+0"

[[deps.Libgpg_error_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "c333716e46366857753e273ce6a69ee0945a6db9"
uuid = "7add5ba3-2f88-524e-9cd5-f83b8a55f7b8"
version = "1.42.0+0"

[[deps.Libiconv_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "f9557a255370125b405568f9767d6d195822a175"
uuid = "94ce4f54-9a6c-5748-9c1c-f9c7231a4531"
version = "1.17.0+0"

[[deps.Libmount_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "dae976433497a2f841baadea93d27e68f1a12a97"
uuid = "4b2f31a3-9ecc-558c-b454-b3730dcb73e9"
version = "2.39.3+0"

[[deps.Libtiff_jll]]
deps = ["Artifacts", "JLLWrappers", "JpegTurbo_jll", "LERC_jll", "Libdl", "XZ_jll", "Zlib_jll", "Zstd_jll"]
git-tree-sha1 = "2da088d113af58221c52828a80378e16be7d037a"
uuid = "89763e89-9b03-5906-acba-b20f662cd828"
version = "4.5.1+1"

[[deps.Libuuid_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "0a04a1318df1bf510beb2562cf90fb0c386f58c4"
uuid = "38a345b3-de98-5d2b-a5d3-14cd9215e700"
version = "2.39.3+1"

[[deps.LinearAlgebra]]
deps = ["Libdl", "OpenBLAS_jll", "libblastrampoline_jll"]
uuid = "37e2e46d-f89d-539d-b4ee-838fcccc9c8e"

[[deps.LogExpFunctions]]
deps = ["DocStringExtensions", "IrrationalConstants", "LinearAlgebra"]
git-tree-sha1 = "18144f3e9cbe9b15b070288eef858f71b291ce37"
uuid = "2ab3a3ac-af41-5b50-aa03-7779005ae688"
version = "0.3.27"

    [deps.LogExpFunctions.extensions]
    LogExpFunctionsChainRulesCoreExt = "ChainRulesCore"
    LogExpFunctionsChangesOfVariablesExt = "ChangesOfVariables"
    LogExpFunctionsInverseFunctionsExt = "InverseFunctions"

    [deps.LogExpFunctions.weakdeps]
    ChainRulesCore = "d360d2e6-b24c-11e9-a2a3-2a2ae2dbcce4"
    ChangesOfVariables = "9e997f8a-9a97-42d5-a9f1-ce6bfc15e2c0"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.Logging]]
uuid = "56ddb016-857b-54e1-b83d-db4d58db5568"

[[deps.LoggingExtras]]
deps = ["Dates", "Logging"]
git-tree-sha1 = "c1dd6d7978c12545b4179fb6153b9250c96b0075"
uuid = "e6f89c97-d47a-5376-807f-9c37f3926c36"
version = "1.0.3"

[[deps.MIMEs]]
git-tree-sha1 = "65f28ad4b594aebe22157d6fac869786a255b7eb"
uuid = "6c6e2e6c-3030-632d-7369-2d6c69616d65"
version = "0.1.4"

[[deps.MacroTools]]
deps = ["Markdown", "Random"]
git-tree-sha1 = "2fa9ee3e63fd3a4f7a9a4f4744a52f4856de82df"
uuid = "1914dd2f-81c6-5fcd-8719-6d5c9610ff09"
version = "0.5.13"

[[deps.Markdown]]
deps = ["Base64"]
uuid = "d6f4376e-aef5-505a-96c1-9c027394607a"

[[deps.MbedTLS]]
deps = ["Dates", "MbedTLS_jll", "MozillaCACerts_jll", "NetworkOptions", "Random", "Sockets"]
git-tree-sha1 = "c067a280ddc25f196b5e7df3877c6b226d390aaf"
uuid = "739be429-bea8-5141-9913-cc70e7f3736d"
version = "1.1.9"

[[deps.MbedTLS_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "c8ffd9c3-330d-5841-b78e-0817d7145fa1"
version = "2.28.2+1"

[[deps.Measures]]
git-tree-sha1 = "c13304c81eec1ed3af7fc20e75fb6b26092a1102"
uuid = "442fdcdd-2543-5da2-b0f3-8c86c306513e"
version = "0.3.2"

[[deps.Missings]]
deps = ["DataAPI"]
git-tree-sha1 = "ec4f7fbeab05d7747bdf98eb74d130a2a2ed298d"
uuid = "e1d29d7a-bbdc-5cf2-9ac0-f12de2c33e28"
version = "1.2.0"

[[deps.Mmap]]
uuid = "a63ad114-7e13-5084-954f-fe012c677804"

[[deps.MozillaCACerts_jll]]
uuid = "14a3606d-f60d-562e-9121-12d972cd8159"
version = "2023.1.10"

[[deps.NaNMath]]
deps = ["OpenLibm_jll"]
git-tree-sha1 = "0877504529a3e5c3343c6f8b4c0381e57e4387e4"
uuid = "77ba4419-2d1f-58cd-9bb1-8ffee604a2e3"
version = "1.0.2"

[[deps.NetworkOptions]]
uuid = "ca575930-c2e3-43a9-ace4-1e988b2c1908"
version = "1.2.0"

[[deps.Ogg_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "887579a3eb005446d514ab7aeac5d1d027658b8f"
uuid = "e7412a2a-1a6e-54c0-be00-318e2571c051"
version = "1.3.5+1"

[[deps.OpenBLAS_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Libdl"]
uuid = "4536629a-c528-5b80-bd46-f80d51c5b363"
version = "0.3.23+2"

[[deps.OpenLibm_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "05823500-19ac-5b8b-9628-191a04bc5112"
version = "0.8.1+2"

[[deps.OpenSSL]]
deps = ["BitFlags", "Dates", "MozillaCACerts_jll", "OpenSSL_jll", "Sockets"]
git-tree-sha1 = "38cb508d080d21dc1128f7fb04f20387ed4c0af4"
uuid = "4d8831e6-92b7-49fb-bdf8-b643e874388c"
version = "1.4.3"

[[deps.OpenSSL_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "3da7367955dcc5c54c1ba4d402ccdc09a1a3e046"
uuid = "458c3c95-2e84-50aa-8efc-19380b2a3a95"
version = "3.0.13+1"

[[deps.Opus_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "51a08fb14ec28da2ec7a927c4337e4332c2a4720"
uuid = "91d4177d-7536-5919-b921-800302f37372"
version = "1.3.2+0"

[[deps.OrderedCollections]]
git-tree-sha1 = "dfdf5519f235516220579f949664f1bf44e741c5"
uuid = "bac558e1-5e72-5ebc-8fee-abe8a469f55d"
version = "1.6.3"

[[deps.PCRE2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "efcefdf7-47ab-520b-bdef-62a2eaa19f15"
version = "10.42.0+1"

[[deps.Parsers]]
deps = ["Dates", "PrecompileTools", "UUIDs"]
git-tree-sha1 = "8489905bcdbcfac64d1daa51ca07c0d8f0283821"
uuid = "69de0a69-1ddd-5017-9359-2bf0b02dc9f0"
version = "2.8.1"

[[deps.Pipe]]
git-tree-sha1 = "6842804e7867b115ca9de748a0cf6b364523c16d"
uuid = "b98c9c47-44ae-5843-9183-064241ee97a0"
version = "1.3.0"

[[deps.Pixman_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "JLLWrappers", "LLVMOpenMP_jll", "Libdl"]
git-tree-sha1 = "64779bc4c9784fee475689a1752ef4d5747c5e87"
uuid = "30392449-352a-5448-841d-b1acce4e97dc"
version = "0.42.2+0"

[[deps.Pkg]]
deps = ["Artifacts", "Dates", "Downloads", "FileWatching", "LibGit2", "Libdl", "Logging", "Markdown", "Printf", "REPL", "Random", "SHA", "Serialization", "TOML", "Tar", "UUIDs", "p7zip_jll"]
uuid = "44cfe95a-1eb2-52ea-b672-e2afdf69b78f"
version = "1.10.0"

[[deps.PlotThemes]]
deps = ["PlotUtils", "Statistics"]
git-tree-sha1 = "1f03a2d339f42dca4a4da149c7e15e9b896ad899"
uuid = "ccf2f8ad-2431-5c83-bf29-c5338b663b6a"
version = "3.1.0"

[[deps.PlotUtils]]
deps = ["ColorSchemes", "Colors", "Dates", "PrecompileTools", "Printf", "Random", "Reexport", "Statistics"]
git-tree-sha1 = "7b1a9df27f072ac4c9c7cbe5efb198489258d1f5"
uuid = "995b91a9-d308-5afd-9ec6-746e21dbc043"
version = "1.4.1"

[[deps.Plots]]
deps = ["Base64", "Contour", "Dates", "Downloads", "FFMPEG", "FixedPointNumbers", "GR", "JLFzf", "JSON", "LaTeXStrings", "Latexify", "LinearAlgebra", "Measures", "NaNMath", "Pkg", "PlotThemes", "PlotUtils", "PrecompileTools", "Printf", "REPL", "Random", "RecipesBase", "RecipesPipeline", "Reexport", "RelocatableFolders", "Requires", "Scratch", "Showoff", "SparseArrays", "Statistics", "StatsBase", "UUIDs", "UnicodeFun", "UnitfulLatexify", "Unzip"]
git-tree-sha1 = "442e1e7ac27dd5ff8825c3fa62fbd1e86397974b"
uuid = "91a5bcdd-55d7-5caf-9e0b-520d859cae80"
version = "1.40.4"

    [deps.Plots.extensions]
    FileIOExt = "FileIO"
    GeometryBasicsExt = "GeometryBasics"
    IJuliaExt = "IJulia"
    ImageInTerminalExt = "ImageInTerminal"
    UnitfulExt = "Unitful"

    [deps.Plots.weakdeps]
    FileIO = "5789e2e9-d7fb-5bc7-8068-2c6fae9b9549"
    GeometryBasics = "5c1252a2-5f33-56bf-86c9-59e7332b4326"
    IJulia = "7073ff75-c697-5162-941a-fcdaad2a7d2a"
    ImageInTerminal = "d8c32880-2388-543b-8c61-d9f865259254"
    Unitful = "1986cc42-f94f-5a68-af5c-568840ba703d"

[[deps.PlutoUI]]
deps = ["AbstractPlutoDingetjes", "Base64", "ColorTypes", "Dates", "FixedPointNumbers", "Hyperscript", "HypertextLiteral", "IOCapture", "InteractiveUtils", "JSON", "Logging", "MIMEs", "Markdown", "Random", "Reexport", "URIs", "UUIDs"]
git-tree-sha1 = "ab55ee1510ad2af0ff674dbcced5e94921f867a9"
uuid = "7f904dfe-b85e-4ff6-b463-dae2292396a8"
version = "0.7.59"

[[deps.PrecompileTools]]
deps = ["Preferences"]
git-tree-sha1 = "5aa36f7049a63a1528fe8f7c3f2113413ffd4e1f"
uuid = "aea7be01-6a6a-4083-8856-8a6e6704d82a"
version = "1.2.1"

[[deps.Preferences]]
deps = ["TOML"]
git-tree-sha1 = "9306f6085165d270f7e3db02af26a400d580f5c6"
uuid = "21216c6a-2e73-6563-6e65-726566657250"
version = "1.4.3"

[[deps.Printf]]
deps = ["Unicode"]
uuid = "de0858da-6303-5e67-8744-51eddeeeb8d7"

[[deps.Qt6Base_jll]]
deps = ["Artifacts", "CompilerSupportLibraries_jll", "Fontconfig_jll", "Glib_jll", "JLLWrappers", "Libdl", "Libglvnd_jll", "OpenSSL_jll", "Vulkan_Loader_jll", "Xorg_libSM_jll", "Xorg_libXext_jll", "Xorg_libXrender_jll", "Xorg_libxcb_jll", "Xorg_xcb_util_cursor_jll", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_keysyms_jll", "Xorg_xcb_util_renderutil_jll", "Xorg_xcb_util_wm_jll", "Zlib_jll", "libinput_jll", "xkbcommon_jll"]
git-tree-sha1 = "37b7bb7aabf9a085e0044307e1717436117f2b3b"
uuid = "c0090381-4147-56d7-9ebc-da0b1113ec56"
version = "6.5.3+1"

[[deps.REPL]]
deps = ["InteractiveUtils", "Markdown", "Sockets", "Unicode"]
uuid = "3fa0cd96-eef1-5676-8a61-b3b8758bbffb"

[[deps.Random]]
deps = ["SHA"]
uuid = "9a3f8284-a2c9-5f02-9a11-845980a1fd5c"

[[deps.RecipesBase]]
deps = ["PrecompileTools"]
git-tree-sha1 = "5c3d09cc4f31f5fc6af001c250bf1278733100ff"
uuid = "3cdcf5f2-1ef4-517c-9805-6587b60abb01"
version = "1.3.4"

[[deps.RecipesPipeline]]
deps = ["Dates", "NaNMath", "PlotUtils", "PrecompileTools", "RecipesBase"]
git-tree-sha1 = "45cf9fd0ca5839d06ef333c8201714e888486342"
uuid = "01d81517-befc-4cb6-b9ec-a95719d0359c"
version = "0.6.12"

[[deps.Reexport]]
git-tree-sha1 = "45e428421666073eab6f2da5c9d310d99bb12f9b"
uuid = "189a3867-3050-52da-a836-e630ba90ab69"
version = "1.2.2"

[[deps.RelocatableFolders]]
deps = ["SHA", "Scratch"]
git-tree-sha1 = "ffdaf70d81cf6ff22c2b6e733c900c3321cab864"
uuid = "05181044-ff0b-4ac5-8273-598c1e38db00"
version = "1.0.1"

[[deps.Requires]]
deps = ["UUIDs"]
git-tree-sha1 = "838a3a4188e2ded87a4f9f184b4b0d78a1e91cb7"
uuid = "ae029012-a4dd-5104-9daa-d747884805df"
version = "1.3.0"

[[deps.SHA]]
uuid = "ea8e919c-243c-51af-8825-aaa63cd721ce"
version = "0.7.0"

[[deps.Scratch]]
deps = ["Dates"]
git-tree-sha1 = "3bac05bc7e74a75fd9cba4295cde4045d9fe2386"
uuid = "6c6a2e73-6563-6170-7368-637461726353"
version = "1.2.1"

[[deps.Serialization]]
uuid = "9e88b42a-f829-5b0c-bbe9-9e923198166b"

[[deps.Showoff]]
deps = ["Dates", "Grisu"]
git-tree-sha1 = "91eddf657aca81df9ae6ceb20b959ae5653ad1de"
uuid = "992d4aef-0814-514b-bc4d-f2e9a6c4116f"
version = "1.0.3"

[[deps.SimpleBufferStream]]
git-tree-sha1 = "874e8867b33a00e784c8a7e4b60afe9e037b74e1"
uuid = "777ac1f9-54b0-4bf8-805c-2214025038e7"
version = "1.1.0"

[[deps.Sockets]]
uuid = "6462fe0b-24de-5631-8697-dd941f90decc"

[[deps.SortingAlgorithms]]
deps = ["DataStructures"]
git-tree-sha1 = "66e0a8e672a0bdfca2c3f5937efb8538b9ddc085"
uuid = "a2af1166-a08f-5f64-846c-94a0d3cef48c"
version = "1.2.1"

[[deps.SparseArrays]]
deps = ["Libdl", "LinearAlgebra", "Random", "Serialization", "SuiteSparse_jll"]
uuid = "2f01184e-e22b-5df5-ae63-d93ebab69eaf"
version = "1.10.0"

[[deps.Statistics]]
deps = ["LinearAlgebra", "SparseArrays"]
uuid = "10745b16-79ce-11e8-11f9-7d13ad32a3b2"
version = "1.10.0"

[[deps.StatsAPI]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1ff449ad350c9c4cbc756624d6f8a8c3ef56d3ed"
uuid = "82ae8749-77ed-4fe6-ae5f-f523153014b0"
version = "1.7.0"

[[deps.StatsBase]]
deps = ["DataAPI", "DataStructures", "LinearAlgebra", "LogExpFunctions", "Missings", "Printf", "Random", "SortingAlgorithms", "SparseArrays", "Statistics", "StatsAPI"]
git-tree-sha1 = "5cf7606d6cef84b543b483848d4ae08ad9832b21"
uuid = "2913bbd2-ae8a-5f71-8c99-4fb6c76f3a91"
version = "0.34.3"

[[deps.SuiteSparse_jll]]
deps = ["Artifacts", "Libdl", "libblastrampoline_jll"]
uuid = "bea87d4a-7f5b-5778-9afe-8cc45184846c"
version = "7.2.1+1"

[[deps.TOML]]
deps = ["Dates"]
uuid = "fa267f1f-6049-4f14-aa54-33bafae1ed76"
version = "1.0.3"

[[deps.Tar]]
deps = ["ArgTools", "SHA"]
uuid = "a4e569a6-e804-4fa4-b0f3-eef7a1d5b13e"
version = "1.10.0"

[[deps.TensorCore]]
deps = ["LinearAlgebra"]
git-tree-sha1 = "1feb45f88d133a655e001435632f019a9a1bcdb6"
uuid = "62fd8b95-f654-4bbd-a8a5-9c27f68ccd50"
version = "0.1.1"

[[deps.Test]]
deps = ["InteractiveUtils", "Logging", "Random", "Serialization"]
uuid = "8dfed614-e22c-5e08-85e1-65c5234f0b40"

[[deps.TranscodingStreams]]
git-tree-sha1 = "71509f04d045ec714c4748c785a59045c3736349"
uuid = "3bb67fe8-82b1-5028-8e26-92a6c54297fa"
version = "0.10.7"
weakdeps = ["Random", "Test"]

    [deps.TranscodingStreams.extensions]
    TestExt = ["Test", "Random"]

[[deps.Tricks]]
git-tree-sha1 = "eae1bb484cd63b36999ee58be2de6c178105112f"
uuid = "410a4b4d-49e4-4fbc-ab6d-cb71b17b3775"
version = "0.1.8"

[[deps.URIs]]
git-tree-sha1 = "67db6cc7b3821e19ebe75791a9dd19c9b1188f2b"
uuid = "5c2747f8-b7ea-4ff2-ba2e-563bfd36b1d4"
version = "1.5.1"

[[deps.UUIDs]]
deps = ["Random", "SHA"]
uuid = "cf7118a7-6976-5b1a-9a39-7adc72f591a4"

[[deps.Unicode]]
uuid = "4ec0a83e-493e-50e2-b9ac-8f72acf5a8f5"

[[deps.UnicodeFun]]
deps = ["REPL"]
git-tree-sha1 = "53915e50200959667e78a92a418594b428dffddf"
uuid = "1cfade01-22cf-5700-b092-accc4b62d6e1"
version = "0.4.1"

[[deps.Unitful]]
deps = ["Dates", "LinearAlgebra", "Random"]
git-tree-sha1 = "3c793be6df9dd77a0cf49d80984ef9ff996948fa"
uuid = "1986cc42-f94f-5a68-af5c-568840ba703d"
version = "1.19.0"

    [deps.Unitful.extensions]
    ConstructionBaseUnitfulExt = "ConstructionBase"
    InverseFunctionsUnitfulExt = "InverseFunctions"

    [deps.Unitful.weakdeps]
    ConstructionBase = "187b0558-2788-49d3-abe0-74a17ed4e7c9"
    InverseFunctions = "3587e190-3f89-42d0-90ee-14403ec27112"

[[deps.UnitfulLatexify]]
deps = ["LaTeXStrings", "Latexify", "Unitful"]
git-tree-sha1 = "e2d817cc500e960fdbafcf988ac8436ba3208bfd"
uuid = "45397f5d-5981-4c77-b2b3-fc36d6e9b728"
version = "1.6.3"

[[deps.Unzip]]
git-tree-sha1 = "ca0969166a028236229f63514992fc073799bb78"
uuid = "41fe7b60-77ed-43a1-b4f0-825fd5a5650d"
version = "0.2.0"

[[deps.Vulkan_Loader_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Wayland_jll", "Xorg_libX11_jll", "Xorg_libXrandr_jll", "xkbcommon_jll"]
git-tree-sha1 = "2f0486047a07670caad3a81a075d2e518acc5c59"
uuid = "a44049a8-05dd-5a78-86c9-5fde0876e88c"
version = "1.3.243+0"

[[deps.Wayland_jll]]
deps = ["Artifacts", "EpollShim_jll", "Expat_jll", "JLLWrappers", "Libdl", "Libffi_jll", "Pkg", "XML2_jll"]
git-tree-sha1 = "7558e29847e99bc3f04d6569e82d0f5c54460703"
uuid = "a2964d1f-97da-50d4-b82a-358c7fce9d89"
version = "1.21.0+1"

[[deps.Wayland_protocols_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "93f43ab61b16ddfb2fd3bb13b3ce241cafb0e6c9"
uuid = "2381bf8a-dfd0-557d-9999-79630e7b1b91"
version = "1.31.0+0"

[[deps.XML2_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libiconv_jll", "Zlib_jll"]
git-tree-sha1 = "532e22cf7be8462035d092ff21fada7527e2c488"
uuid = "02c8fc9c-b97f-50b9-bbe4-9be30ff0a78a"
version = "2.12.6+0"

[[deps.XSLT_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Libgcrypt_jll", "Libgpg_error_jll", "Libiconv_jll", "Pkg", "XML2_jll", "Zlib_jll"]
git-tree-sha1 = "91844873c4085240b95e795f692c4cec4d805f8a"
uuid = "aed1982a-8fda-507f-9586-7b0439959a61"
version = "1.1.34+0"

[[deps.XZ_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "ac88fb95ae6447c8dda6a5503f3bafd496ae8632"
uuid = "ffd25f8a-64ca-5728-b0f7-c24cf3aae800"
version = "5.4.6+0"

[[deps.Xorg_libICE_jll]]
deps = ["Libdl", "Pkg"]
git-tree-sha1 = "e5becd4411063bdcac16be8b66fc2f9f6f1e8fe5"
uuid = "f67eecfb-183a-506d-b269-f58e52b52d7c"
version = "1.0.10+1"

[[deps.Xorg_libSM_jll]]
deps = ["Libdl", "Pkg", "Xorg_libICE_jll"]
git-tree-sha1 = "4a9d9e4c180e1e8119b5ffc224a7b59d3a7f7e18"
uuid = "c834827a-8449-5923-a945-d239c165b7dd"
version = "1.2.3+0"

[[deps.Xorg_libX11_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxcb_jll", "Xorg_xtrans_jll"]
git-tree-sha1 = "afead5aba5aa507ad5a3bf01f58f82c8d1403495"
uuid = "4f6342f7-b3d2-589e-9d20-edeb45f2b2bc"
version = "1.8.6+0"

[[deps.Xorg_libXau_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "6035850dcc70518ca32f012e46015b9beeda49d8"
uuid = "0c0b7dd1-d40b-584c-a123-a41640f87eec"
version = "1.0.11+0"

[[deps.Xorg_libXcursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXfixes_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "12e0eb3bc634fa2080c1c37fccf56f7c22989afd"
uuid = "935fb764-8cf2-53bf-bb30-45bb1f8bf724"
version = "1.2.0+4"

[[deps.Xorg_libXdmcp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "34d526d318358a859d7de23da945578e8e8727b7"
uuid = "a3789734-cfe1-5b06-b2d0-1dd0d9d62d05"
version = "1.1.4+0"

[[deps.Xorg_libXext_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "b7c0aa8c376b31e4852b360222848637f481f8c3"
uuid = "1082639a-0dae-5f34-9b06-72781eeb8cb3"
version = "1.3.4+4"

[[deps.Xorg_libXfixes_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "0e0dc7431e7a0587559f9294aeec269471c991a4"
uuid = "d091e8ba-531a-589c-9de9-94069b037ed8"
version = "5.0.3+4"

[[deps.Xorg_libXi_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXfixes_jll"]
git-tree-sha1 = "89b52bc2160aadc84d707093930ef0bffa641246"
uuid = "a51aa0fd-4e3c-5386-b890-e753decda492"
version = "1.7.10+4"

[[deps.Xorg_libXinerama_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll"]
git-tree-sha1 = "26be8b1c342929259317d8b9f7b53bf2bb73b123"
uuid = "d1454406-59df-5ea1-beac-c340f2130bc3"
version = "1.1.4+4"

[[deps.Xorg_libXrandr_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libXext_jll", "Xorg_libXrender_jll"]
git-tree-sha1 = "34cea83cb726fb58f325887bf0612c6b3fb17631"
uuid = "ec84b674-ba8e-5d96-8ba1-2a689ba10484"
version = "1.5.2+4"

[[deps.Xorg_libXrender_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libX11_jll"]
git-tree-sha1 = "19560f30fd49f4d4efbe7002a1037f8c43d43b96"
uuid = "ea2f1a96-1ddc-540d-b46f-429655e07cfa"
version = "0.9.10+4"

[[deps.Xorg_libpthread_stubs_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "8fdda4c692503d44d04a0603d9ac0982054635f9"
uuid = "14d82f49-176c-5ed1-bb49-ad3f5cbd8c74"
version = "0.1.1+0"

[[deps.Xorg_libxcb_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "XSLT_jll", "Xorg_libXau_jll", "Xorg_libXdmcp_jll", "Xorg_libpthread_stubs_jll"]
git-tree-sha1 = "b4bfde5d5b652e22b9c790ad00af08b6d042b97d"
uuid = "c7cfdc94-dc32-55de-ac96-5a1b8d977c5b"
version = "1.15.0+0"

[[deps.Xorg_libxkbfile_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libX11_jll"]
git-tree-sha1 = "730eeca102434283c50ccf7d1ecdadf521a765a4"
uuid = "cc61e674-0454-545c-8b26-ed2c68acab7a"
version = "1.1.2+0"

[[deps.Xorg_xcb_util_cursor_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xcb_util_image_jll", "Xorg_xcb_util_jll", "Xorg_xcb_util_renderutil_jll"]
git-tree-sha1 = "04341cb870f29dcd5e39055f895c39d016e18ccd"
uuid = "e920d4aa-a673-5f3a-b3d7-f755a4d47c43"
version = "0.1.4+0"

[[deps.Xorg_xcb_util_image_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "0fab0a40349ba1cba2c1da699243396ff8e94b97"
uuid = "12413925-8142-5f55-bb0e-6d7ca50bb09b"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_libxcb_jll"]
git-tree-sha1 = "e7fd7b2881fa2eaa72717420894d3938177862d1"
uuid = "2def613f-5ad1-5310-b15b-b15d46f528f5"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_keysyms_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "d1151e2c45a544f32441a567d1690e701ec89b00"
uuid = "975044d2-76e6-5fbe-bf08-97ce7c6574c7"
version = "0.4.0+1"

[[deps.Xorg_xcb_util_renderutil_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "dfd7a8f38d4613b6a575253b3174dd991ca6183e"
uuid = "0d47668e-0667-5a69-a72c-f761630bfb7e"
version = "0.3.9+1"

[[deps.Xorg_xcb_util_wm_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Xorg_xcb_util_jll"]
git-tree-sha1 = "e78d10aab01a4a154142c5006ed44fd9e8e31b67"
uuid = "c22f9ab0-d5fe-5066-847c-f4bb1cd4e361"
version = "0.4.1+1"

[[deps.Xorg_xkbcomp_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_libxkbfile_jll"]
git-tree-sha1 = "330f955bc41bb8f5270a369c473fc4a5a4e4d3cb"
uuid = "35661453-b289-5fab-8a00-3d9160c6a3a4"
version = "1.4.6+0"

[[deps.Xorg_xkeyboard_config_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Xorg_xkbcomp_jll"]
git-tree-sha1 = "691634e5453ad362044e2ad653e79f3ee3bb98c3"
uuid = "33bec58e-1273-512f-9401-5d533626f822"
version = "2.39.0+0"

[[deps.Xorg_xtrans_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e92a1a012a10506618f10b7047e478403a046c77"
uuid = "c5fb5394-a638-5e4d-96e5-b29de1b5cf10"
version = "1.5.0+0"

[[deps.Zlib_jll]]
deps = ["Libdl"]
uuid = "83775a58-1f1d-513f-b197-d71354ab007a"
version = "1.2.13+1"

[[deps.Zstd_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "e678132f07ddb5bfa46857f0d7620fb9be675d3b"
uuid = "3161d3a3-bdf6-5164-811a-617609db77b4"
version = "1.5.6+0"

[[deps.eudev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "gperf_jll"]
git-tree-sha1 = "431b678a28ebb559d224c0b6b6d01afce87c51ba"
uuid = "35ca27e7-8b34-5b7f-bca9-bdc33f59eb06"
version = "3.2.9+0"

[[deps.fzf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl"]
git-tree-sha1 = "a68c9655fbe6dfcab3d972808f1aafec151ce3f8"
uuid = "214eeab7-80f7-51ab-84ad-2988db7cef09"
version = "0.43.0+0"

[[deps.gperf_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3516a5630f741c9eecb3720b1ec9d8edc3ecc033"
uuid = "1a1c6b14-54f6-533d-8383-74cd7377aa70"
version = "3.1.1+0"

[[deps.libaom_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "3a2ea60308f0996d26f1e5354e10c24e9ef905d4"
uuid = "a4ae2306-e953-59d6-aa16-d00cac43593b"
version = "3.4.0+0"

[[deps.libass_jll]]
deps = ["Artifacts", "Bzip2_jll", "FreeType2_jll", "FriBidi_jll", "HarfBuzz_jll", "JLLWrappers", "Libdl", "Pkg", "Zlib_jll"]
git-tree-sha1 = "5982a94fcba20f02f42ace44b9894ee2b140fe47"
uuid = "0ac62f75-1d6f-5e53-bd7c-93b484bb37c0"
version = "0.15.1+0"

[[deps.libblastrampoline_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850b90-86db-534c-a0d3-1478176c7d93"
version = "5.8.0+1"

[[deps.libevdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "141fe65dc3efabb0b1d5ba74e91f6ad26f84cc22"
uuid = "2db6ffa8-e38f-5e21-84af-90c45d0032cc"
version = "1.11.0+0"

[[deps.libfdk_aac_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "daacc84a041563f965be61859a36e17c4e4fcd55"
uuid = "f638f0a6-7fb0-5443-88ba-1cc74229b280"
version = "2.0.2+0"

[[deps.libinput_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "eudev_jll", "libevdev_jll", "mtdev_jll"]
git-tree-sha1 = "ad50e5b90f222cfe78aa3d5183a20a12de1322ce"
uuid = "36db933b-70db-51c0-b978-0f229ee0e533"
version = "1.18.0+0"

[[deps.libpng_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Zlib_jll"]
git-tree-sha1 = "d7015d2e18a5fd9a4f47de711837e980519781a4"
uuid = "b53b4c65-9356-5827-b1ea-8c7a1a84506f"
version = "1.6.43+1"

[[deps.libvorbis_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Ogg_jll", "Pkg"]
git-tree-sha1 = "b910cb81ef3fe6e78bf6acee440bda86fd6ae00c"
uuid = "f27f6e37-5d2b-51aa-960f-b287f2bc3b7a"
version = "1.3.7+1"

[[deps.mtdev_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "814e154bdb7be91d78b6802843f76b6ece642f11"
uuid = "009596ad-96f7-51b1-9f1b-5ce2d5e8a71e"
version = "1.1.6+0"

[[deps.nghttp2_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "8e850ede-7688-5339-a07c-302acd2aaf8d"
version = "1.52.0+1"

[[deps.p7zip_jll]]
deps = ["Artifacts", "Libdl"]
uuid = "3f19e933-33d8-53b3-aaab-bd5110c3b7a0"
version = "17.4.0+2"

[[deps.x264_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "4fea590b89e6ec504593146bf8b988b2c00922b2"
uuid = "1270edf5-f2f9-52d2-97e9-ab00b5d0237a"
version = "2021.5.5+0"

[[deps.x265_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg"]
git-tree-sha1 = "ee567a171cce03570d77ad3a43e90218e38937a9"
uuid = "dfaa095f-4041-5dcd-9319-2fabd8486b76"
version = "3.5.0+0"

[[deps.xkbcommon_jll]]
deps = ["Artifacts", "JLLWrappers", "Libdl", "Pkg", "Wayland_jll", "Wayland_protocols_jll", "Xorg_libxcb_jll", "Xorg_xkeyboard_config_jll"]
git-tree-sha1 = "9c304562909ab2bab0262639bd4f444d7bc2be37"
uuid = "d8fb68d0-12a3-5cfd-a85a-d49703b185fd"
version = "1.4.1+1"
"""

# ╔═╡ Cell order:
# ╟─e052fe10-c756-4128-b5e4-7c19d2d389a6
# ╠═fd9d5608-8d6f-4904-af06-c93d5b9a05bd
# ╠═04db8de6-4267-4df9-b38a-217a46cce879
# ╠═44993ac9-2ed3-4f05-8fdf-245cee136561
# ╟─c38eae7c-b9f5-4591-83ce-4e1f6f04c8ed
# ╠═a5c8b870-fa82-4e5c-9934-625b755e95c8
# ╠═07c26752-1429-4cd7-bc71-f4fb3fde7b1e
# ╟─712f4850-77cc-4e5e-b74f-bfa075e337fe
# ╠═147aec25-366c-4450-a0b4-088f1fe4ead9
# ╟─7dd51cc1-d3aa-4d94-b1d8-01205eb1aa23
# ╟─9b7096ea-1e09-4bf3-9e0f-a34561f7ccbf
# ╟─84814e56-cc4a-4f20-ab92-72100d5bf0d0
# ╠═2bc32094-7427-4e02-a76d-fc77722fc280
# ╠═76cf69bb-9af7-431c-98fe-e427d40f6552
# ╟─afb89567-1d3f-4a39-bfa8-d99cffeac6dd
# ╠═50f84de9-9092-43f9-b0c2-abe76e40342f
# ╟─c285ca6d-8d40-41cd-b643-be8704798336
# ╟─61e8e74c-ccb6-46dc-a5da-46f188784f70
# ╟─b2cba02d-08da-445b-ad34-8971bfcdd882
# ╠═86bf6b83-d989-47e7-a758-64b458d826a8
# ╠═18254000-21c0-407d-8d31-080321b51c6d
# ╟─3915c82f-f36b-4ca7-9357-c799f25bc0a1
# ╠═22643b1c-8b92-4e46-b071-7f4c99019949
# ╟─2b296a5c-6a47-459c-8b56-8100850c59ca
# ╠═903a72db-3357-4018-ba28-2dfb5905714b
# ╟─f3216d7f-6978-4123-8e72-f28b859dea26
# ╠═fd606282-1ff2-4d78-a120-2ac660b219fa
# ╠═0345834c-1818-4023-9de8-04f732269dba
# ╠═07692bc0-e289-4ac3-af4e-2ee7ec29b6a2
# ╠═861c4fb6-13cf-43db-b8f2-4cff0142dcd2
# ╟─fa897a98-3304-417f-9074-7c84a1c20356
# ╟─de757199-f956-4e69-bdc0-7910b971ef34
# ╠═212b15dd-5ba1-49f8-8fbc-c28e8dcfa69a
# ╟─82cee6e0-36f5-4ff4-9517-3c4ce95c4fa5
# ╟─43ff924f-8b74-4947-ae4d-5b6c28b58eba
# ╟─8b5a4da4-85c7-41f7-8568-cd89d6ab2eea
# ╠═0ba4cf06-cb37-44a9-8963-92fa3c660dc4
# ╠═2a5ea354-efa2-45e5-9f60-9c4e3d77540a
# ╟─2271c992-8f4a-4ad8-a04a-111404ed2576
# ╠═597c6890-f1d3-4f2d-9031-9c834b7d9968
# ╟─6a21606e-8d0d-4502-97c4-9895b39f9d21
# ╠═3f613597-abc4-4b38-ae03-238856fc89a3
# ╟─499dcecc-8c23-4570-83b5-b776517f665a
# ╟─1f5b5feb-6ccd-4e6c-8979-1c688409924d
# ╠═7429fc5e-7275-4ae1-876b-0fa83a5714c0
# ╟─04f92c6d-b580-44ea-a9c8-41e13fe02d25
# ╟─00000000-0000-0000-0000-000000000001
# ╟─00000000-0000-0000-0000-000000000002
