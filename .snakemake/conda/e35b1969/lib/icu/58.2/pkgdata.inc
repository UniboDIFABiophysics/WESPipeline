GENCCODE_ASSEMBLY_TYPE=-a gcc
SO=so
SOBJ=so
A=a
LIBPREFIX=lib
LIB_EXT_ORDER=.58.2
COMPILE=gcc -I/home/PERSONALE/eugenio.fonzi2/WESPipeline/.snakemake/conda/e35b1969/include -D_REENTRANT  -DU_HAVE_ELF_H=1 -DU_HAVE_ATOMIC=1 -DU_HAVE_STRTOD_L=1  -DU_ATTRIBUTE_DEPRECATED= -m64 -std=c99 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -c
LIBFLAGS=-I/home/PERSONALE/eugenio.fonzi2/WESPipeline/.snakemake/conda/e35b1969/include -DPIC -fPIC
GENLIB=gcc -m64 -std=c99 -Wall -pedantic -Wshadow -Wpointer-arith -Wmissing-prototypes -Wwrite-strings   -Wl,-rpath,/home/PERSONALE/eugenio.fonzi2/WESPipeline/.snakemake/conda/e35b1969/lib -L/home/PERSONALE/eugenio.fonzi2/WESPipeline/.snakemake/conda/e35b1969/lib  -shared -Wl,-Bsymbolic
LDICUDTFLAGS=-nodefaultlibs -nostdlib
LD_SONAME=-Wl,-soname -Wl,
RPATH_FLAGS=
BIR_LDFLAGS=-Wl,-Bsymbolic
AR=ar
ARFLAGS=r
RANLIB=ranlib
INSTALL_CMD=/usr/bin/install -c
