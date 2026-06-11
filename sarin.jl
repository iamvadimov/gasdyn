using Random
using Printf

# Структура для представлення Лагранжевої частинки
mutable struct GasParticle
    x::Float64
    y::Float64
    z::Float64
    mass::Float64      # Маса зарину в частинці (кг)
    is_trapped::Bool   # Прапорець перебування у підвалі (застійна зона)
end

# Фізичні, географічні та хімічні константи
const g = 9.81              # Прискорення вільного падіння (м/с²)
const ρ_air = 1.204         # Густина повітря (кг/м³)
const ρ_sarin_vapor = 5.85  # Густина пари зарину при стандартних умовах (кг/м³)
const P_sat = 387.0         # Тиск насиченої пари зарину при 20°C (Па)
const M_sarin = 0.1401      # Молярна маса зарину (кг/моль)
const R = 8.314             # Універсальна газова стала

# Параметри кінетики розпаду зарину
const E_a = 41500.0          # Енергія активації гідролізу (Дж/моль)
const A_h = 120.0            # Масштабний коефіцієнт гідролізу
const α_uv = 1.5e-6          # Коефіцієнт фотодеструкції (м²/Вт·с)

# Конфігурація геометрії сцени
const BASEMENT_Z = -3.0     # Глибина підвалу (м)
const BASEMENT_X_MIN = 4.0;  const BASEMENT_X_MAX = 8.0
const BASEMENT_Y_MIN = 4.0;  const BASEMENT_Y_MAX = 8.0

const STACK_X_MIN = 5.0;     const STACK_X_MAX = 7.0
const STACK_Y_MIN = 5.0;     const STACK_Y_MAX = 7.0
const STACK_H = 15.0        # Висота будівлі з ефектом димової труби

"""
Функція обчислення сумарної константи деструкції зарину (1/с)
"""
function calculate_destruction_rate(T_k::Float64, RH::Float64, I_solar::Float64)
    k_hydro = A_h * (RH / 100.0) * exp(-E_a / (R * T_k))
    k_photo = α_uv * I_solar
    return k_hydro + k_photo
end

"""
Функція обчислення швидкості випаровування зарину (кг/с)
"""
function calculate_evaporation(A::Float64, U_wind::Float64, T_k::Float64)
    K_m = 5.5e-6 + 4.2e-6 * U_wind
    dm_dt = K_m * A * (P_sat * M_sarin / (R * T_k))
    return max(0.0, dm_dt)
end

"""
Основна функция кроку симуляції Лагранжевих частинок з урахуванням деструкції маси
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
    # Швидкість гравітаційного затікання важкого газу
    w_g = -sqrt(2.0 * g * 0.5 * (ρ_sarin_vapor - ρ_air) / ρ_air)
    
    # Виштовхувальна швидкість ефекту димової труби
    w_stack = 0.6 * sqrt(2.0 * g * STACK_H * max(0.0, T_in - T_out) / T_in)

    # Видалення дегрованих частинок із немізерною масою
    filter!(p -> p.mass > 1e-8, particles)

    for p in particles
        # Перевірка геометрії підвалу
        in_basement_geo = (p.x >= BASEMENT_X_MIN && p.x <= BASEMENT_X_MAX &&
                           p.y >= BASEMENT_Y_MIN && p.y <= BASEMENT_Y_MAX &&
                           p.z <= 0.0)

        if in_basement_geo
            p.is_trapped = true
        end

        # Зменшення маси частинки внаслідок деструкції
        current_k = p.is_trapped ? k_decay_basement : k_decay
        p.mass *= exp(-current_k * dt)

        if p.is_trapped
            # Дифузія та гравітація всередині підвалу (ефект кишені)
            p.x += σ_turb * sqrt(dt) * randn()
            p.y += σ_turb * sqrt(dt) * randn()
            p.z += w_g * dt + σ_turb * sqrt(dt) * randn()
            p.z = max(BASEMENT_Z, p.z)
            
            # Вентиляційний вихід частинки з підвалу
            if rand() < (1.0 - exp(-λ_v * dt))
                p.is_trapped = false
                p.z = 0.1
            end
        else
            # Рух у міському каньйоні з логарифмічним профілем вітру
            wind_profile = p.z > 0.0 ? U_mean * log(max(1.0, p.z) / 0.1) / log(10.0 / 0.1) : 0.1 * U_mean
            
            p.x += wind_profile * dt + σ_turb * sqrt(dt) * randn()
            p.y += σ_turb * sqrt(dt) * randn()
            
            dz_eff = w_g * dt + σ_turb * sqrt(dt) * randn()
            
            # Перевірка потрапляння в зону інфільтрації будівлі
            in_stack_zone = (p.x >= STACK_X_MIN && p.x <= STACK_X_MAX &&
                             p.y >= STACK_Y_MIN && p.y <= STACK_Y_MAX &&
                             p.z >= 0.0 && p.z <= STACK_H)
            if in_stack_zone
                dz_eff += w_stack * dt
            end
            p.z += dz_eff

            # Ідеальне відбиття від підстильної поверхні
            if p.z < 0.0 && !in_basement_geo
                p.z = -p.z
            end
        end
    end
end

"""
Демонстраційний сценарій локального витоку зарину
"""
function run_advanced_simulation()
    println("--- Старт симуляції експрес-оцінки з урахуванням деструкції ---")
    
    # Метеорологічні параметри середовища
    T_ambient = 298.15  # 25 °C на вулиці
    RH_outdoor = 50.0   # 50% відносна вологість ззовні
    I_solar = 800.0     # Інтенсивна інсоляція (ясний день, Вт/м²)
    
    T_basement = 288.15 # 15 °C в підвальному приміщенні
    RH_basement = 90.0  # 90% відносна вологість у підвалі
    
    # Попередній розрахунок кінетичних констант
    k_outdoor = calculate_destruction_rate(T_ambient, RH_outdoor, I_solar)
    k_indoor = calculate_destruction_rate(T_basement, RH_basement, 0.0) # Сонце відсутнє
    
    @printf("Розрахункова константа деструкції (Вулиця): %.6f c⁻¹\n", k_outdoor)
    @printf("Розрахункова константа деструкції (Підвал): %.6f c⁻¹\n", k_indoor)

    # Параметри джерела: пролив ~3 літрів рідкої фази (маса ~2.97 кг, площа ~2 м²)
    M_total = 2.97 
    A_spill = 2.0
    particles = GasParticle[]
    
    dt = 0.5     # Крок інтегрування за часом (с)
    t_end = 60.0 # Загальний час аналізу (с)
    time_elapsed = 0.0
    
    σ_turb = 0.25 # Інтенсивність мікротурбулентних пульсацій
    λ_v = 0.02   # Кратність повітрообміну підвалу
    
    while time_elapsed < t_end
        # Динамічне випаровування рідкої фази за крок dt
        dm = calculate_evaporation(A_spill, 3.0, T_ambient) * dt
        if M_total > 0
            actual_dm = min(dm, M_total)
            M_total -= actual_dm
            
            # Дискретизація маси пари на ансамбль Лагранжевих частинок
            num_new = ceil(Int, actual_dm * 100) 
            for _ in 1:num_new
                push!(particles, GasParticle(2.0, 5.0, 0.0, actual_dm/num_new, false))
            end
        end
        
        # Обчислення нового стану системи частинок
        update_particles_with_decay!(particles, dt, 3.0, σ_turb, λ_v, 295.15, T_ambient, k_outdoor, k_indoor)
        
        # Інтегральний аналіз маси
        total_mass_active = sum(p -> p.mass, particles, init=0.0)
        mass_in_basement = sum(p -> p.is_trapped ? p.mass : 0.0, particles, init=0.0)
        
        if rem(time_elapsed, 10.0) == 0
            @printf("Час: %4.1f с | Активна маса в повітрі: %5.3f кг | Маса у підвалах: %5.3f кг | Часток: %d\n", 
                    time_elapsed, total_mass_active, mass_in_basement, length(particles))
        end
        time_elapsed += dt
    end
    println("--- Розрахунок сценарію завершено успішно ---")
end

run_advanced_simulation()
