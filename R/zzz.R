
# Load package-specific options

.onLoad <- function(libname, pkgname) {
    .set_xmcutil_opts(pkgname)
}
