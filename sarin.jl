using Random
using Printf

# Structure representing a Lagrangian particle
mutable struct GasParticle
    x::Float64
    y::Float64
    z::Float64
    mass::Float64      # Sarin mass in the particle (kg)
    is_trapped::Bool   # Flag indicating presence in a basement (stagnation zone)
end

# Physical, geographical, and chemical constants
const g = 9.81              # Gravitational acceleration (m/s²)
const ρ_air = 1.204         # Air density (kg/m³)
const ρ_sarin_vapor = 5.85  # Sarin vapor density under standard conditions (kg/m³)
const P_sat = 387.0         # Saturated vapor pressure of Sarin at 20°C (Pa)
const M_sarin = 0.1401      # Molar mass of Sarin (kg/mol)
const R = 8.314             # Universal gas constant

# Sarin degradation kinetics parameters
const E_a = 41500.0          # Hydrolysis activation energy (J/mol)
const A_h = 120.0            # Hydrolysis scaling coefficient
const α_uv = 1.5e-6          # Photodegradation coefficient (m²/W·s)

# Scene geometry configuration
const BASEMENT_Z = -3.0     # Basement depth (m)
const BASEMENT_X_MIN = 4.0;  const BASEMENT_X_MAX = 8.0
const BASEMENT_Y_MIN = 4.0;  const BASEMENT_Y_MAX = 8.0

const STACK_X_MIN = 5.0;     const STACK_X_MAX = 7.0
const STACK_Y_MIN = 5.0;     const STACK_Y_MAX = 7.0
const STACK_H = 15.0        # Building height with stack effect

"""
Function for calculating the total Sarin destruction rate constant (1/s)
"""
function calculate_destruction_rate(T_k::Float64, RH::Float64, I_solar::Float64)
    k_hydro = A_h * (RH / 100.0) * exp(-E_a / (R * T_k))
    k_photo = α_uv * I_solar
    return k_hydro + k_photo
end

"""
Function for calculating the Sarin evaporation rate (kg/s)
"""
function calculate_evaporation(A::Float64, U_wind::Float64, T_k::Float64)
    K_m = 5.5e-6 + 4.2e-6 * U_wind
    dm_dt = K_m * A * (P_sat * M_sarin / (R * T_k))
    return max(0.0, dm_dt)
end

"""
Main Lagrangian particle simulation step with mass degradation
"""
function update_particles_with_decay!(
    particles::Vector{GasParticle}, 
    dt::Float64, 
    U_mean::Float64, 
    σ_turb::Float64, 
    λ_v::Float64, 
    T_in::Float64, 
    T_out::Float64,
    k_decay::Float64,
    k_decay_basement::Float64
)
    # Gravity-driven settling velocity of heavy gas
    w_g = -sqrt(2.0 * g * 0.5 * (ρ_sarin_vapor - ρ_air) / ρ_air)
    
    # Upward velocity induced by the stack effect
    w_stack = 0.6 * sqrt(2.0 * g * STACK_H * max(0.0, T_in - T_out) / T_in)

    # Remove degraded particles with negligible mass
    filter!(p -> p.mass > 1e-8, particles)

    for p in particles
        # Check basement geometry
        in_basement_geo = (p.x >= BASEMENT_X_MIN && p.x <= BASEMENT_X_MAX &&
                           p.y >= BASEMENT_Y_MIN && p.y <= BASEMENT_Y_MAX &&
                           p.z <= 0.0)

        if in_basement_geo
            p.is_trapped = true
        end

        # Reduce particle mass due to degradation
        current_k = p.is_trapped ? k_decay_basement : k_decay
        p.mass *= exp(-current_k * dt)

        if p.is_trapped
            # Diffusion and gravity inside the basement (pocket effect)
            p.x += σ_turb * sqrt(dt) * randn()
            p.y += σ_turb * sqrt(dt) * randn()
            p.z += w_g * dt + σ_turb * sqrt(dt) * randn()
            p.z = max(BASEMENT_Z, p.z)
            
            # Ventilation-driven particle escape from basement
            if rand() < (1.0 - exp(-λ_v * dt))
                p.is_trapped = false
                p.z = 0.1
            end
        else
            # Motion in the urban canyon with logarithmic wind profile
            wind_profile = p.z > 0.0 ? U_mean * log(max(1.0, p.z) / 0.1) / log(10.0 / 0.1) : 0.1 * U_mean
            
            p.x += wind_profile * dt + σ_turb * sqrt(dt) * randn()
            p.y += σ_turb * sqrt(dt) * randn()
            
            dz_eff = w_g * dt + σ_turb * sqrt(dt) * randn()
            
            # Check entry into the building infiltration zone
            in_stack_zone = (p.x >= STACK_X_MIN && p.x <= STACK_X_MAX &&
                             p.y >= STACK_Y_MIN && p.y <= STACK_Y_MAX &&
                             p.z >= 0.0 && p.z <= STACK_H)
            if in_stack_zone
                dz_eff += w_stack * dt
            end
            p.z += dz_eff

            # Perfect reflection from the ground surface
            if p.z < 0.0 && !in_basement_geo
                p.z = -p.z
            end
        end
    end
end

"""
Demonstration scenario of a local Sarin release
"""
function run_advanced_simulation()
    println("--- Start of rapid assessment simulation with degradation ---")
    
    # Environmental meteorological parameters
    T_ambient = 298.15  # 25 °C outdoors
    RH_outdoor = 50.0   # 50% outdoor relative humidity
    I_solar = 800.0     # High solar irradiance (clear day, W/m²)
    
    T_basement = 288.15 # 15 °C inside the basement
    RH_basement = 90.0  # 90% relative humidity in the basement
    
    # Preliminary calculation of kinetic constants
    k_outdoor = calculate_destruction_rate(T_ambient, RH_outdoor, I_solar)
    k_indoor = calculate_destruction_rate(T_basement, RH_basement, 0.0) # No solar radiation
    
    @printf("Calculated destruction constant (Outdoor): %.6f c⁻¹\n", k_outdoor)
    @printf("Calculated destruction constant (Basement): %.6f c⁻¹\n", k_indoor)

    # Source parameters: spill of ~3 liters of liquid phase (mass ~2.97 kg, area ~2 m²)
    M_total = 2.97 
    A_spill = 2.0
    particles = GasParticle[]
    
    dt = 0.5     # Time integration step (s)
    t_end = 60.0 # Total simulation time (s)
    time_elapsed = 0.0
    
    σ_turb = 0.25 # Intensity of micro-turbulent fluctuations
    λ_v = 0.02   # Basement air exchange rate
    
    while time_elapsed < t_end
        # Dynamic evaporation of the liquid phase during time step dt
        dm = calculate_evaporation(A_spill, 3.0, T_ambient) * dt
        if M_total > 0
            actual_dm = min(dm, M_total)
            M_total -= actual_dm
            
            # Discretization of vapor mass into an ensemble of Lagrangian particles
            num_new = ceil(Int, actual_dm * 100) 
            for _ in 1:num_new
                push!(particles, GasParticle(2.0, 5.0, 0.0, actual_dm/num_new, false))
            end
        end
        
        # Compute the new state of the particle system
        update_particles_with_decay!(particles, dt, 3.0, σ_turb, λ_v, 295.15, T_ambient, k_outdoor, k_indoor)
        
        # Integral mass analysis
        total_mass_active = sum(p -> p.mass, particles, init=0.0)
        mass_in_basement = sum(p -> p.is_trapped ? p.mass : 0.0, particles, init=0.0)
        
        if rem(time_elapsed, 10.0) == 0
            @printf("Time: %4.1f с | Active mass in air: %5.3f кг | Mass in basements: %5.3f кг | Particles: %d\n", 
                    time_elapsed, total_mass_active, mass_in_basement, length(particles))
        end
        time_elapsed += dt
    end
    println("--- Scenario simulation completed successfully ---")
end

run_advanced_simulation()
