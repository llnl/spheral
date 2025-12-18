-- -*- lua -*-
help([[
Loads Spheral/YYYY.MM.p release
 Commands: spheral
 Modifies: Add /usr/gapps/Spheral/$SYS_TYPE/YYYY.MM.p to the path in order to use the YYYY.MM.p release version of the Spheral modeling tools.
]])
whatis("simulation")
whatis("Spheral release YYYY.MM.p")
whatis([[ Loads Spheral/YYYY.MM.p
 Commands: spheral
 Modifies: Add /usr/gapps/Spheral/$SYS_TYPE/YYYY.MM.p to the path in order to use the YYYY.MM.p release version of the Spheral modeling tools.
]])

prepend_path("PATH", "/usr/gapps/Spheral/" .. os.getenv("SYS_TYPE") .. "/vYYYY.MM.p/bin")
