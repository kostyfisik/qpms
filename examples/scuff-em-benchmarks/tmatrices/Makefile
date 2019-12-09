tmatrices: tmatrices_veryrough tmatrices_rough tmatrices_normal tmatrices_fine

omegalists_scuff: omegalist2_scuff omegalist1_scuff

omegalist1_scuff: omegalist1_eV	
	../../../misc/omega_eV2scuff.py -o omegalist1_scuff omegalist1_eV 

omegalist2_scuff: omegalist2_eV	
	../../../misc/omega_eV2scuff.py -o omegalist2_scuff omegalist2_eV 

tmatrices_veryrough: cylinderAu_r100_h50_lMax3_veryrough.TMatrix cylinderAg_r30_h30_lMax3_veryrough.TMatrix

tmatrices_rough: cylinderAu_r100_h50_lMax3_rough.TMatrix cylinderAg_r30_h30_lMax3_rough.TMatrix

tmatrices_fine: cylinderAu_r100_h50_lMax3_fine.TMatrix cylinderAg_r30_h30_lMax3_fine.TMatrix

tmatrices_normal: cylinderAu_r100_h50_lMax3.TMatrix cylinderAg_r30_h30_lMax3.TMatrix

cylinderAu_r100_h50.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r100_h50.msh/' -e 's/__THEMATERIAL__/LDGold/' cylinder.scuffgeo.template > $@

cylinderAu_r100_h50_veryrough.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r100_h50_veryrough.msh/' -e 's/__THEMATERIAL__/LDGold/' cylinder.scuffgeo.template > $@

cylinderAu_r100_h50_rough.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r100_h50_rough.msh/' -e 's/__THEMATERIAL__/LDGold/' cylinder.scuffgeo.template > $@

cylinderAu_r100_h50_fine.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r100_h50_fine.msh/' -e 's/__THEMATERIAL__/LDGold/' cylinder.scuffgeo.template > $@

cylinderAg_r30_h30.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r30_h30.msh/' -e 's/__THEMATERIAL__/LDSilver/' cylinder.scuffgeo.template > $@

cylinderAg_r30_h30_veryrough.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r30_h30_veryrough.msh/' -e 's/__THEMATERIAL__/LDSilver/' cylinder.scuffgeo.template > $@

cylinderAg_r30_h30_rough.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r30_h30_rough.msh/' -e 's/__THEMATERIAL__/LDSilver/' cylinder.scuffgeo.template > $@

cylinderAg_r30_h30_fine.scuffgeo: cylinder.scuffgeo.template
	sed -e 's/__THEMESHFILE__/cylinder_r30_h30_fine.msh/' -e 's/__THEMATERIAL__/LDSilver/' cylinder.scuffgeo.template > $@

cylinderAu_r100_h50_lMax3.TMatrix: omegalist1_scuff cylinder_r100_h50.msh cylinderAu_r100_h50.scuffgeo
	scuff-tmatrix --geometry cylinderAu_r100_h50.scuffgeo --OmegaFile omegalist1_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAu_r100_h50_lMax3_veryrough.TMatrix: omegalist1_scuff cylinder_r100_h50_veryrough.msh cylinderAu_r100_h50_veryrough.scuffgeo
	scuff-tmatrix --geometry cylinderAu_r100_h50_veryrough.scuffgeo --OmegaFile omegalist1_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAu_r100_h50_lMax3_rough.TMatrix: omegalist1_scuff cylinder_r100_h50_rough.msh cylinderAu_r100_h50_rough.scuffgeo
	scuff-tmatrix --geometry cylinderAu_r100_h50_rough.scuffgeo --OmegaFile omegalist1_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAu_r100_h50_lMax3_fine.TMatrix: omegalist1_scuff cylinder_r100_h50_fine.msh cylinderAu_r100_h50_fine.scuffgeo
	scuff-tmatrix --geometry cylinderAu_r100_h50_fine.scuffgeo --OmegaFile omegalist1_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAg_r30_h30_lMax3_veryrough.TMatrix: omegalist2_scuff cylinder_r30_h30_veryrough.msh cylinderAg_r30_h30_veryrough.scuffgeo
	scuff-tmatrix --geometry cylinderAg_r30_h30_veryrough.scuffgeo --OmegaFile omegalist2_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAg_r30_h30_lMax3_rough.TMatrix: omegalist2_scuff cylinder_r30_h30_rough.msh cylinderAg_r30_h30_rough.scuffgeo
	scuff-tmatrix --geometry cylinderAg_r30_h30_rough.scuffgeo --OmegaFile omegalist2_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAg_r30_h30_lMax3_fine.TMatrix: omegalist2_scuff cylinder_r30_h30_fine.msh cylinderAg_r30_h30_fine.scuffgeo
	scuff-tmatrix --geometry cylinderAg_r30_h30_fine.scuffgeo --OmegaFile omegalist2_scuff --FileBase $(@:.TMatrix=) --LMax 3

cylinderAg_r30_h30_lMax3.TMatrix: omegalist2_scuff cylinder_r30_h30.msh cylinderAg_r30_h30.scuffgeo
	scuff-tmatrix --geometry cylinderAg_r30_h30.scuffgeo --OmegaFile omegalist2_scuff --FileBase $(@:.TMatrix=) --LMax 3
