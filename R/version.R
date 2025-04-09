print(paste("R", getRversion()))
print("-------------")
for (package_name in sort(loadedNamespaces())) {
  print(paste(package_name, packageVersion(package_name)))
}
