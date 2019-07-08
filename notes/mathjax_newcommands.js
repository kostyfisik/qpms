MathJax.Hub.Config({
    TeX: {
        Macros: {
            vect: ["{\\mathbf{#1}}",1],
            abs: ["{\\left|{#1}\\right|}",1],
            ud: "{\\mathrm{d}}",
	    pr: ["{\\left({#1}\\right)}", 1], // parentheses to save typing
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
	    wfmtreg: "{\\widetilde{\\vect{M}}^{(1)}}", // regular magnetic wave
	    wfetreg: "{\\widetilde{\\vect{N}}^{(1)}}", // regular electric wave
	    wfmtout: "{\\widetilde{\\vect{M}}^{(3)}}", // outgoing magnetic wave
	    wfetout: "{\\widetilde{\\vect{N}}^{(3)}}", // outgoing electric wave
	    wcmtreg: "{q}", // regular magnetic wave coeff
	    wcetreg: "{p}", // regular electric wave coeff
	    wcmtout: "{b}", // outgoing magnetic wave coeff
	    wcetout: "{a}" // outgoing electric wave coeff
        }
    }
});
