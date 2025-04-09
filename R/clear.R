rm(list = ls())

# This does not work too well consider restarting R environment if package clashes.
#
# basicPackages <-  c("compiler","graphics","tools","rstudioapi","utils","grDevices","stats","datasets","methods","base")
# pkgs <- setdiff(loadedNamespaces(), basicPackages)
# 
# 
# while (length(pkgs)>0){
#   pkgs <- setdiff(loadedNamespaces(), c("compiler","graphics","tools","rstudioapi","utils","grDevices","stats","datasets","methods","base"))
#   for (pkg in rev(pkgs)) {  # Reverse order to handle dependencies first
#     try({
#       unloadNamespace(pkg)  # Unload namespace
#       detach(paste0("package:", pkg), character.only = TRUE, unload = TRUE)
#     }, silent = TRUE)
#   }
#   setdiff(loadedNamespaces(),basicPackages)
# }

#setdiff(loadedNamespaces(),basicPackages)
#rm(list = ls())
