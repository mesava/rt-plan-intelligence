# Планируемая модель данных исследования

Все ключи обезличенные. Исходные пациентские DICOM/XCC/PDF не публикуются.

- **План:** anonymous_plan_key, rtplan_uid_hash, technique, site, TPS/version, linac_model, MLC_model, energy, fractions, metrics_version.
- **QA:** measured_gamma_percent, gamma_DD_percent, gamma_DTA_mm, gamma_global_local, gamma_low_dose_cutoff, detector, QA software/version, measurement geometry, calibration version, normalization, reviewer_status.
- **Луч/CP:** beam index, MU, gantry angle, MLC positions, jaws, gap/aperture/travel metrics.
- **Сводный анализ:** группировка по технике, машине, локализации и протоколу измерения.

Требования: не допускать утечки данных одного пациента между training/test; контролировать версии алгоритмов; не смешивать разные gamma-протоколы без учёта их различий. Модель не обучена и не интегрирована в v0.1.0.
