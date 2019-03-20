enable_language (Fortran)
include(FortranCInterface)

FortranCInterface_HEADER(amos_mangling.h
	MACRO_NAMESPACE "AMOS_"
	SYMBOL_NAMESPACE "amos_"
	SYMBOLS zbesj zbesy zbesh zbesi zbesk zsqrt
	)

add_library(amos 
	dgamln.f zabs.f zasyi.f zbinu.f zdiv.f zmlri.f zshch.f zunhj.f zunk1.f
	d1mach.f zacai.f zbesh.f zbknu.f zexp.f zmlt.f zsqrt.f zunik.f zunk2.f
	i1mach.f zacon.f zbesj.f zbuni.f zkscl.f zrati.f zs1s2.f zuni1.f zuoik.f
	xerror.f zairy.f zbesy.f zbunk.f zlog.f zseri.f zuchk.f zuni2.f zwrsk.f
	)

target_include_directories (amos PUBLIC	${CMAKE_CURRENT_SOURCE_DIR}) 	

