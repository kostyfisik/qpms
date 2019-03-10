#!/usr/bin/env python3
from qpms.symmetries import point_group_info

codestring = "#include <qpms/groups.h>\n"
for name, info in point_group_info.items():
    codestring += 'const qpms_finite_group_t QPMS_FINITE_GROUP_%s = ' %name
    codestring += info.generate_c_source()
    codestring += ";\n\n"

print(codestring)


