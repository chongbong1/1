# Анализ FMM для 2D логарифмического потенциала

## Проблема

Текущая реализация имеет **КРИТИЧЕСКУЮ ОШИБКУ** в формуле M2L (Multipole-to-Local translation).

## Правильные формулы для 2D логарифмического ядра

### 1. Фундаментальное разложение

Для логарифмического потенциала в 2D:
```
φ(z) = -log|z - z₀|
```

Ключевое разложение в комплексных переменных:
```
-log|z - z₀| = -log|z| - Re[Σ_{k=1}^∞ (z₀/z)^k / k]  для |z₀| < |z|
```

### 2. P2M: Multipole expansion (правильно!)

Мультипольные моменты относительно центра z_c:
```
M₀ = Σᵢ qᵢ                                    (монополь)
M_k = -Σᵢ qᵢ · (zᵢ/rscale)^k / k              (k ≥ 1)
```

✓ **Это реализовано ПРАВИЛЬНО**

### 3. M2L: Translation (ОШИБКА НАЙДЕНА!)

#### Текущая НЕПРАВИЛЬНАЯ формула:

```fortran
! НЕПРАВИЛЬНО:
L₀ = M₀ * log|z₀/rscale| + Σ_{k=1}^p M_k

L_j = [-M₀/j + Σ_{k=1}^p M_k * C(j+k-1,k-1)] * (1/z₀)^j
```

#### ПРАВИЛЬНАЯ формула из литературы:

Для перевода мультипольного разложения из источника в локальное разложение в цели:

**Вектор смещения**: z₀ = z_target - z_source (от источника к цели)

**Формулы M2L:**

```
L₀ = M₀ * (-log|z₀|) + Σ_{k=1}^∞ M_k * (-1/z₀)^(-k)
   = M₀ * (-log|z₀|) + Σ_{k=1}^∞ M_k * z₀^k

L_j = -M₀/(j·z₀^j) + Σ_{k=1}^∞ M_k · C(j+k-1, k-1) / z₀^(j+k)     (j ≥ 1)
```

где C(n,k) - биномиальные коэффициенты.

### 4. Ключевые ошибки в текущем коде:

#### Ошибка #1: Неправильный знак для L₀

```fortran
! НЕПРАВИЛЬНО:
sys%cells(ti,tj)%L(0) = sys%cells(ti,tj)%L(0) + &
                       sys%cells(si,sj)%M(0) * log(abs(z0_scaled))

! ПРАВИЛЬНО (должен быть МИНУС):
sys%cells(ti,tj)%L(0) = sys%cells(ti,tj)%L(0) - &
                       sys%cells(si,sj)%M(0) * log(abs(z0_scaled))
```

#### Ошибка #2: Неправильная формула для вклада M_k в L₀

```fortran
! НЕПРАВИЛЬНО:
do k = 1, sys%p_max
  sys%cells(ti,tj)%L(0) = sys%cells(ti,tj)%L(0) + sys%cells(si,sj)%M(k)
end do

! ПРАВИЛЬНО (должен быть z₀^k, не просто M_k):
z0_pow = z0_scaled
do k = 1, sys%p_max
  sys%cells(ti,tj)%L(0) = sys%cells(ti,tj)%L(0) + sys%cells(si,sj)%M(k) * z0_pow
  z0_pow = z0_pow * z0_scaled
end do
```

#### Ошибка #3: Неправильный знак степени в L_j

```fortran
! НЕПРАВИЛЬНО (используется 1/z₀^j):
z0_pow = z0_inv_scaled  ! (1/z0)
do j = 1, sys%p_max
  temp_sum = -sys%cells(si,sj)%M(0) / real(j, dp)
  do k = 1, sys%p_max
    temp_sum = temp_sum + sys%cells(si,sj)%M(k) * binom_coef
  end do
  sys%cells(ti,tj)%L(j) = sys%cells(ti,tj)%L(j) + temp_sum * z0_pow
  z0_pow = z0_pow * z0_inv_scaled
end do

! ПРАВИЛЬНО (нужно 1/z₀^(j+k)):
z0_inv_pow_j = z0_inv_scaled  ! (1/z0)^1
do j = 1, sys%p_max
  ! Монопольный вклад: -M₀/(j·z₀^j)
  temp_sum = -sys%cells(si,sj)%M(0) / real(j, dp) * z0_inv_pow_j

  ! Вклад высших моментов: Σ M_k · C(j+k-1,k-1) / z₀^(j+k)
  z0_inv_pow_jk = z0_inv_pow_j * z0_inv_scaled  ! (1/z0)^(j+1)
  do k = 1, sys%p_max
    binom_coef = sys%binomial(j+k-1, k-1)
    temp_sum = temp_sum + sys%cells(si,sj)%M(k) * binom_coef * z0_inv_pow_jk
    z0_inv_pow_jk = z0_inv_pow_jk * z0_inv_scaled  ! (1/z0)^(j+k)
  end do

  sys%cells(ti,tj)%L(j) = sys%cells(ti,tj)%L(j) + temp_sum
  z0_inv_pow_j = z0_inv_pow_j * z0_inv_scaled
end do
```

## Источники правильных формул

1. **Greengard & Rokhlin (1987)** - оригинальная статья FMM
2. **fmm2d library (Flatiron Institute)** - современная эталонная реализация
3. **Beatson & Greengard FMM course** - учебные материалы

## Проверка формул

### Физический смысл:

- Потенциал логарифмический: φ = -q·log(r)
- Сила: F = -∇φ = q·r/r²
- M2L должен корректно переводить разложение в дальней зоне

### Математическая консистентность:

Правильные формулы получаются из разложения:
```
φ(r_target) = -Σᵢ qᵢ log|r_target - r_source,i|
            = -Σᵢ qᵢ log|z₀ - zᵢ|     где z₀ = r_target - r_source_center
```

Используя разложение log в комплексной плоскости и группируя члены по степеням z₀.

## Рекомендации

1. **Исправить M2L формулы** согласно анализу выше
2. **Увеличить p_order до 30-40** для высокой точности (у тебя сейчас 20)
3. **Проверить направление z₀** - должно быть от источника к цели
4. **Добавить диагностику**: вывод значений L_k после M2L для проверки

## Ожидаемая точность

После исправления при p_order = 30-40:
- Относительная ошибка силы: < 10⁻¹² (машинная точность)
- Для p_order = 20: ошибка ~ 10⁻⁸ - 10⁻¹⁰
