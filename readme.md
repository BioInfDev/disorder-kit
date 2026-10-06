Для анализа степени неупорядоченности белковых последовательностей были выбраны следующие метрики: AIUPred , pLDDT и counterharmonic mean.
AIUPred использует AI модель для поиска участков неупорядоченности в белковой последовательности. Для каждой аминокислоты рассчитывается score - вероятность того, что конкретный аминокислотный остаток принадлежит к внутренне неупорядоченному региону.
pLDDT - мера локальной уверенности модели в локальных межатомных расстояниях вокруг остатка (по Calpha). Метрика используется для оценки качества предсказанной тритичной структуры белка. Однако была показана умеренная корреляция pLDDT и IDRs (doi.org/10.1093/nar/gkad928): "we show the correlation of the disorder content, i.e. the fraction of disordered residues in the protein sequence, between DisProt and AlphaFold2 when different pLDDT thresholds are selected. DisProt and AlphaFold2 correlate well when the pLDDT threshold is between 70 and 90 with a maximum correlation at pLDDT = 80 (Pearson's correlation 0.42, P-value 6.26e−101)". Наличие корреляция позволяет использовать pLDDT для скрининга белковых последовательностей на предмет наличия IDRs. 

AIUPred и pLDDT - per-residue метрики, поэтому для характеристики неупорядоченности белка в целом рассчитывается процент аминокислотных остатков со значениями AIUPred score > 0.5 и pLDDT < 80 соответственно.

Counterharmonic mean - контргармоническое среднее (взвешенное среднее арифметическое (wi = xi)) - позволяет обобщить AIUPred и pLDDT. Контргармоническое срденее учитывает вклад величин пропорционально их значению. 

