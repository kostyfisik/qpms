/*! \file scatsystem.h
 * \brief Modern interface for finite lattice calculations, including symmetries.
 */
#ifndef QPMS_SCATSYSTEM_H
#define QPMS_SCATSYSTEM_H

/// A T-matrix.
/** In the future, I might rather use a more abstract approach in which T-matrix
 *  is a mapping (function) of the field expansion coefficients.
 *  So the interface might change.
 *  For now, let me stick to the square dense matrix representation.
 */
typedef struct qpms_tmatrix_t {
	const qpms_vswf_set_spec_t *spec; ///< NOT owned by qpms_tmatrix_t by default.
	complex double *m; ///< Matrix elements 
} qpms_tmatrix_t;
/// Returns a pointer to the beginning of the T-matrix row \a rowno.
static inline complex double *qpms_tmatrix_row(qpms_tmatrix_t *t, size_t rowno){
	return t->m + (t->spec->n * rowno);
}
/// NOT IMPLEMENTED Initialises a zero T-matrix.
qpms_tmatrix_t *qpms_tmatrix_init(const qpms_vswf_set_spec_t *bspec);
/// NOT IMPLEMENTED Destroys a T-matrix.
void qpms_tmatrix_free(qpms_tmatrix_t *t);


/// A particle, defined by its T-matrix and position.
typedef struct qpms_particle_t {
	// Does it make sense to ever use other than cartesian coords for this?
	cart3_t pos; ///< Particle position in cartesian coordinates.
	const qpms_tmatrix_t *tmatrix; ///< T-matrix; not owned by qpms_particle_t.
} qpms_particle_t;

#endif //QPMS_SCATSYSTEM_H
