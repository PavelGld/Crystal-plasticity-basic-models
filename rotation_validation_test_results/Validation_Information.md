# Russian

В данной папке приводится информация по теме валидации модели. Валидация прогнозирования кривых напряжение-деформация производилась по статье Shveykin A., Romanov K., Trusov P. Some issues with statistical crystal plasticity models: Description of the effects triggered in FCC crystals by loading with strain-path changes //materials. – 2022. – Т. 15. – №. 19. – С. 6586.

Результаты расчета представленной модели содержатся в файлах Stress_strain_hard_loading.txt и Stress_strain_unaxial_tension.txt. Их сопоставление с графиками из статьи приводится в Stress_strain_verification_curves.xlsx. Наблюдаемое отклонение не превышает 4-5 МПа, что может быть вызвано обычным отклонением от шага интегрирования, либо различий в представительных объемах.

Модель ротации Тейлора проверялась по полюсным фигурам нагружений:
1. Квазиодноосное сжатие
2. Квазиодноосное растяжение
3. Простой сдвиг
4. Плоское сжатие

Полюсные фигуры и .txt кватернионов зерен находятся в соответствующих файлах ("orientations_<deformation_type>_<percent_of_def>pct", "PF_<deformation_type>_<percent_of_def>pct"). Сами полюсные фигуры расчитывались с применением открытого программного обеспечения MTEX (matlab). Для перевода кватернионов в функцию распределения ориентаций и последующего построения полюсных фигур был написать скрипт plot_polefigures.m
Полученные фигуры сопоставлялись с фигурами из статей:

Anand L. Single-crystal elasto-viscoplasticity: application to texture evolution in polycrystalline metals at large strains //Computer methods in applied mechanics and engineering. – 2004. – Т. 193. – №. 48-51. – С. 5359-5383.
Bronkhorst C. A., Kalidindi S. R., Anand L. Polycrystalline plasticity and the evolution of crystallographic texture in FCC metals //Philosophical Transactions of the Royal Society of London. Series A: Physical and Engineering Sciences. – 1992. – Т. 341. – №. 1662. – С. 443-477.

Важно заметить, что хоть в двух вышеупомянутых работах использовалась медь как материал, в верификации рассматривался алюминий, что тем не менее не сильно влияет на итоговые полюсные фигуры, т.к. оба материала имеют ГЦК решетку.

Сопоставление показало хорошее соответствие полюсных фигур с наблюдаемыми экспериментально и численно в других работах, а также сходимость по кривым напряжение-деформация, что доказывает работоспособность представленной в репозитории реализации базовой модели.


# English
This folder contains information on model validation. The prediction of stress-strain curves was validated using the article by Shveykin A., Romanov K., Trusov P. "Some issues with statistical crystal plasticity models: Description of the effects triggered in FCC crystals by loading with strain-path changes" [in Russian]. Materials. – 2022. – Vol. 15. – No. 19. – P. 6586.

The calculation results for the presented model are contained in the files Stress_strain_hard_loading.txt and Stress_strain_unaxial_tension.txt. A comparison of these results with the graphs from the article is provided in Stress_strain_verification_curves.xlsx. The observed deviation does not exceed 4-5 MPa, which may be due to a normal deviation from the integration step or differences in the representative volumes.

The Taylor rotation model was verified using pole figures for the following loading conditions:
1. Quasi-uniaxial compression
2. Quasi-uniaxial tension
3. Simple shear
4. Plane compression

The pole figures and .txt grain quaternions are located in the corresponding files ("orientations_<deformation_type>_<percent_of_def>pct", "PF_<deformation_type>_<percent_of_def>pct"). The pole figures themselves were calculated using the open-source MTEX ​​software (matlab). To convert the quaternions into orientation distribution functions and subsequently plot the pole figures, the plot_polefigures.m script was written.
The resulting figures were compared with those from the following articles:

Anand L. Single-crystal elasto-viscoplasticity: application to texture evolution in polycrystalline metals at large strains // Computer methods in applied mechanics and engineering. – 2004. – Vol. 193. – No. 48-51. – Pp. 5359-5383.
Bronkhorst C. A., Kalidindi S. R., Anand L. Polycrystalline plasticity and the evolution of crystallographic texture in FCC metals // Philosophical Transactions of the Royal Society of London. Series A: Physical and Engineering Sciences. – 1992. – Vol. 341. – No. 1662. – Pp. 443-477.

It is important to note that although the two aforementioned studies used copper as the material, aluminum was considered in the verification study. This, however, does not significantly affect the resulting pole figures, since both materials have an fcc lattice.

The comparison showed good agreement between the pole figures and those observed experimentally and numerically in other studies, as well as convergence in the stress-strain curves, which proves the viability of the implementation of the basic model presented in the repository.