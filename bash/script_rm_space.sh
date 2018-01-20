 #!/*bin/*bash
 # A bash script which removes white characters at the end of each line from all files listed in a given directory. A name of a directory in which files will be modified should be passed as an argument, e.g. bash script_rm_char Documents. If the name of a directory is not listed the script will remove white characters from files listed in 'git status'. Script processes only text files (not binary files and not directories). Script processes files that have a space in the name.
 
 echo "White characters at the end of lines in all files " 
 if [ -n "$1" ]
 then
        echo "in '$1' directory are going to be removed" 
 # Looking for files in given directory"
	MY_FILES=$(ls $1)
	cd $1

 else
 	echo "given by 'git status' are going to be removed"
 # Looking for files printed in 'git status'"
 	MY_FILES=$(git status | \
 	sed -i -n '/^$/,/^$/p' | \ 
 	sed -i -e 's/^.*\(modified:\|new file:\|deleted.*$\)//g' | \ 
	sed -i '/^$/d') 
 fi
 
 if [ -n "$MY_FILES" ] # True if the length of "$MY_FILES" is non-zero. Alternatively the -n could be ignored: if [ "$MY_FILES" ]
 then
	 echo "The edited files are:"
     		for l in $MY_FILES
     		do
         		if file -b "$l" | grep -q 'text'; then
            		sed -i 's/[[:blank:]]*$//'  "$l"
            		echo "$l"
        		fi
     		done 
 else echo "No files to edit"
 fi
 echo "Done"
