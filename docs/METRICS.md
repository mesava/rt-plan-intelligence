# Метрики v0.1.0: определения и ограничения

**Только плановые данные DICOM RTPLAN MLCX.** Все расстояния в миллиметрах на уровне изоцентра, площади в мм². Эти показатели не измеряют gamma и не доказывают доставляемость.

| Поле | Определение |
| --- | --- |
| `control_points` | Число CP на луч/арку |
| `monitor_units` | BeamMeterset из FractionGroupSequence, если задан |
| `aperture_area_mean_mm2` | MU-взвешенное по ΔCMW среднее площади апертуры после clipping по jaws; берётся конечный CP интервала |
| `aperture_area_cv` | SD площадей CP / mean площадей CP, без MU-веса |
| `leaf_gap_mean_mm`, `leaf_gap_median_mm` | Среднее и медиана положительных jaw-clipped зазоров открытых пар по CP |
| `sas5_mu_weighted`, `sas10_mu_weighted` | Доля открытых и экспонированных Y jaws пар со значением 0 < gap < 5 / 10 мм, взвешенная по ΔCMW конечного CP |
| `leaf_travel_total_mm` | Сумма модулей смещений всех лепестков обоих банков между соседними CP (в том числе вне jaws) |
| `leaf_travel_mean_per_transition_mm`, `leaf_travel_max_per_transition_mm` | Среднее и максимум суммарного движения за переход |
| `mean_cp_aperture_area_change_fraction` | Среднее |Aᵢ₊₁−Aᵢ|/(Aᵢ₊₁+Aᵢ) по переходам, где сумма площадей > 0 |

Расчёт площади пары MLCX:

`gap_i = max(0, min(right_i, xjaw_max) - max(left_i, xjaw_min))`

`height_i = max(0, min(y_upper_i, yjaw_max) - max(y_lower_i, yjaw_min))`

`area = Σ gap_i × height_i`

При отсутствии jaws по X принимается неограниченная X-область; по Y — полная высота MLC.

**Важные замечания:**
- SAS по конечному CP — это **приближение** динамической доставки, а не непрерывный интеграл и не обязательно SAS по литературной формуле.
- Total Leaf Travel учитывает лепестки вне jaws, поэтому сравнение с другими реализациями требует совпадения определения.
- `aperture_area_cv` **не является литературным AAV/MCS**. `mean_cp_aperture_area_change_fraction` **не является валидированным Temporal Aperture Variability**.
- Пустые/закрытые поля и нулевая ΔCMW могут давать `null`, а не ноль.
- Для расчёта **скоростей и ускорений** нужна достоверная временная шкала, которой одних CP недостаточно.
- MCS, LSV, AAV, XCC, RTDOSE/RTSTRUCT, gamma prediction и clinical QA verdict **не реализованы** в v0.1.0.
