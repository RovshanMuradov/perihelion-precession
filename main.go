package main

import (
	"fmt"
	"math"
)

// Физические константы
const (
	G  = 6.67430e-11    // Гравитационная постоянная (м³/(кг·с²))
	c  = 299792458      // Скорость света (м/с)
	AU = 1.495978707e11 // Астрономическая единица (м)

	// Параметры Солнца
	M_sun = 1.98892e30 // Масса Солнца (кг)

	// Параметры орбиты Меркурия
	a_mercury = 0.387098 * AU // Большая полуось (м)
	e_mercury = 0.20563       // Эксцентриситет
	T_mercury = 87.969        // Период обращения (дни)
)

// Структура для хранения результатов
type PrecessionResult struct {
	PerOrbit         float64 // Прецессия за один оборот (радианы)
	PerOrbitArcSec   float64 // Прецессия за один оборот (угловые секунды)
	PerCentury       float64 // Прецессия за столетие (угловые секунды)
	OrbitsPerCentury float64 // Количество оборотов за столетие
}

// Вычисление прецессии перигелия по формуле ОТО
func calculateGRPrecession(mass, semiMajorAxis, eccentricity, orbitalPeriodDays float64) PrecessionResult {
	// Формула ОТО для прецессии за один оборот:
	// Δφ = 6πGM / (c²a(1-e²))

	numerator := 6 * math.Pi * G * mass
	denominator := c * c * semiMajorAxis * (1 - eccentricity*eccentricity)

	precessionPerOrbit := numerator / denominator

	// Перевод из радиан в угловые секунды
	// 1 радиан = 206264.806 угловых секунд
	precessionPerOrbitArcSec := precessionPerOrbit * 206264.806

	// Количество оборотов за 100 лет
	daysPerCentury := 365.25 * 100
	orbitsPerCentury := daysPerCentury / orbitalPeriodDays

	// Прецессия за столетие
	precessionPerCentury := precessionPerOrbitArcSec * orbitsPerCentury

	return PrecessionResult{
		PerOrbit:         precessionPerOrbit,
		PerOrbitArcSec:   precessionPerOrbitArcSec,
		PerCentury:       precessionPerCentury,
		OrbitsPerCentury: orbitsPerCentury,
	}
}

// Вычисление вклада от других планет (приближенно)
func calculatePlanetaryPerturbations() float64 {
	// Известные вклады от других планет (угловые секунды за столетие)
	venus := 277.8
	earth := 90.0
	jupiter := 153.6
	otherPlanets := 10.0 // Марс, Сатурн и другие

	return venus + earth + jupiter + otherPlanets
}

func main() {
	fmt.Println("=== Расчет прецессии перигелия Меркурия ===")
	fmt.Println()

	// Вычисляем прецессию по ОТО
	result := calculateGRPrecession(M_sun, a_mercury, e_mercury, T_mercury)

	fmt.Printf("Параметры орбиты Меркурия:\n")
	fmt.Printf("- Большая полуось: %.3f АЕ (%.3e м)\n", a_mercury/AU, a_mercury)
	fmt.Printf("- Эксцентриситет: %.5f\n", e_mercury)
	fmt.Printf("- Период обращения: %.3f дней\n", T_mercury)
	fmt.Printf("- Оборотов за столетие: %.1f\n", result.OrbitsPerCentury)
	fmt.Println()

	fmt.Printf("Прецессия по Общей Теории Относительности:\n")
	fmt.Printf("- За один оборот: %.6e радиан\n", result.PerOrbit)
	fmt.Printf("- За один оборот: %.3f угловых секунд\n", result.PerOrbitArcSec)
	fmt.Printf("- За столетие: %.1f угловых секунд\n", result.PerCentury)
	fmt.Println()

	// Вклад от других планет
	planetaryEffect := calculatePlanetaryPerturbations()
	fmt.Printf("Вклад от других планет (классическая механика): %.1f\"\n", planetaryEffect)
	fmt.Println()

	// Общая предсказанная прецессия
	totalPredicted := planetaryEffect + result.PerCentury
	observedPrecession := 574.1 // Наблюдаемое значение

	fmt.Printf("Сравнение с наблюдениями:\n")
	fmt.Printf("- Наблюдаемая прецессия: %.1f\" за столетие\n", observedPrecession)
	fmt.Printf("- Предсказание классической механики: %.1f\"\n", planetaryEffect)
	fmt.Printf("- Дополнительный вклад ОТО: %.1f\"\n", result.PerCentury)
	fmt.Printf("- Полное предсказание (классика + ОТО): %.1f\"\n", totalPredicted)
	fmt.Printf("- Разница с наблюдениями: %.1f\"\n", math.Abs(totalPredicted-observedPrecession))
	fmt.Println()

	// Проверяем точность
	if math.Abs(result.PerCentury-43.0) < 1.0 {
		fmt.Println("✓ Расчет успешен! Получены знаменитые ~43 угловые секунды!")
	}

	// Дополнительная информация
	fmt.Println("\n=== Интересные факты ===")
	fmt.Printf("- Эффект ОТО составляет всего %.1f%% от общей прецессии\n",
		result.PerCentury/observedPrecession*100)
	fmt.Printf("- За один оборот Меркурия перигелий смещается всего на %.3f угловых секунд\n",
		result.PerOrbitArcSec)
	fmt.Printf("- Это примерно 1/1000 углового размера Луны!\n")
}

// Дополнительные функции для расширения проекта

// Расчет для произвольной планеты
func calculateForPlanet(name string, mass, a, e, T float64) {
	fmt.Printf("\n=== Расчет для %s ===\n", name)
	result := calculateGRPrecession(mass, a, e, T)
	fmt.Printf("Прецессия ОТО за столетие: %.2f угловых секунд\n", result.PerCentury)
}

// Визуализация прецессии (псевдокод)
func visualizePrecession() {
	// Здесь можно добавить код для создания визуализации
	// например, используя библиотеку для построения графиков
	fmt.Println("\n[Здесь могла бы быть визуализация орбиты с прецессией]")
}
