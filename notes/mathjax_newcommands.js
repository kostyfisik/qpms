MathJax.Hub.Config({
    TeX: {
        Macros: {

		//    Abs: ['\\left\\lvert #2 \\right\\rvert_{\\text{#1}}', 2, ""] // optional arg. example
		//    from https://stackoverflow.com/questions/24628668/how-to-define-custom-macros-in-mathjax
            vect: ["{\\mathbf{#1}}",1],
            abs: ["{\\left|{#1}\\right|}",1],
            ud: "{\\mathrm{d}}",
	    pr: ["{\\left({#1}\\right)}", 1], // parentheses to save typing
	    uvec: ["{\\mathbf{\\hat{#1}}}", 1],

	    vsh: "{\\mathbf{A}}", // vector spherical harmonic, general
	    vshD: "\\mathbf{A}^\\dagger", // dual vector spherical harmonic, general
	    vshrad: "{\\mathbf{A}_3}", // vector spherical harmonic radial, general
	    vshrot: "{\\mathbf{A}_1}", // vector spherical harmonic "rotational", general
	    vshgrad: "{\\mathbf{A}_2}", // vector spherical harmonic "gradiental", general
	    vshradD: "{\\mathbf{A}_3}^\\dagger}", // dual vector spherical harmonic radial, general
	    vshrotD: "{\\mathbf{A}_1^\\dagger}", // dual vector spherical harmonic "rotational", general
	    vshgradD: "{\\mathbf{A}_2^\\dagger}", // dual vector spherical harmonic "gradiental", general
	    wfe: "{\\mathbf{N}}", // Electric wave general
	    wfm: "{\\mathbf{M}}", // Magnetic wave general
	    sphbes: "{z}", // General spherical Bessel fun
	    rawLeg: ["{P_{#1}^{#2}}", 2], // "Canonical" associated Legendre polynomial

            // Kristensson's VSWFs, complex version (2014 notes)
            wfkcreg: "{\\vect{v}}", // regular wave
            wfkcout: "{\\vect{u}}", // outgoing wave
            wckcreg: "{a}", // regular wave coeff
            wckcout: "{f}", // outgoing wave coeff
            
	    // Kristensson's VSWFs, real version (2014 book)
            wfkrreg: "{\\vect{v}}", // regular wave
            wfkrout: "{\\vect{u}}", // outgoing wave
            wckrreg: "{a}", // regular wave coeff
            wckrout: "{f}", // outgoing wave coeff

	    // Taylor's VSWFs
	    wfmt: "{\\widetilde{\\vect{M}}}",
	    wfet: "{\\widetilde{\\vect{N}}}",
	    wfmtreg: "{\\widetilde{\\vect{M}}^{(1)}}", // regular magnetic wave
	    wfetreg: "{\\widetilde{\\vect{N}}^{(1)}}", // regular electric wave
	    wfmtout: "{\\widetilde{\\vect{M}}^{(3)}}", // outgoing magnetic wave
	    wfetout: "{\\widetilde{\\vect{N}}^{(3)}}", // outgoing electric wave
	    wcmtreg: "{q}", // regular magnetic wave coeff
	    wcetreg: "{p}", // regular electric wave coeff
	    wcmtout: "{b}", // outgoing magnetic wave coeff
	    wcetout: "{a}", // outgoing electric wave coeff

	    // Reid's VSWFs
	    wfr: "\\mathbf{\\mathcal{W}}",
            wfrreg: "\\mathbf{\\mathcal{W}}^{\\mathrm{reg}}", // regular wave
            wfrout: "\\mathbf{\\mathcal{W}}^{\\mathrm{out}}", // outgoing wave
            wcrreg: "C^\\mathrm{inc}", // regular wave coeff
            wcrout: "C^\\mathrm{scat}", // outgoing wave coeff
        }
    }
});

