# TUTORIAL I: Basic UNIX Commands  
Mon Oct 29 11:11:35 EDT 2018  
By Deborah Velez-Irizarry  
  
This tutorial is a review of some basic UNIX commands such as creating  
directories, saving and read files from the command line and some  
useful commands to keep in mind. These are core skills you  
will use every time you work on HPC.
  
### Connect to HPC system  
There are various ways to connect to HPC.  
   
> Web-based remote desktop environment: [Web site access to HPCC](https://wiki.hpcc.msu.edu/display/ITH/Web+Site+Access+to+HPCC)  
  
> Connect to HPCC by SSH client is avalable at [iCER](https://wiki.hpcc.msu.edu/display/ITH/Connect+to+HPC+System)  
  
Take some time to learn more about MSU HPCC. Some good resources to look over are:  
  
> [Why use HPCC](https://wiki.hpcc.msu.edu/display/ITH/Why+Use+HPCC)  
> [HPC Glossary](https://wiki.hpcc.msu.edu/display/ITH/HPC+Glossary)  
> [HPC entire layout at iCER](https://wiki.hpcc.msu.edu/display/ITH/HPC%27s+entire+layout+at+iCER)  
  
  
### Lets get started.  
  
Today we will work on our first tutorial of UNIX. You will learn how to   
create directories, save and read files from the command line and some  
usefull commands to keep in mind. These are core skills you will use every  
time you work on HPC.  
  
  
### Home Directory  
  
Lets start by looking at your home directory on HPC. When you log in to HPC  
you are automatically taken to your home directory. Your home environment is  
saved under the tilde `~` symbol. This is convenient because you can go to  
your home directory by just using the change directory `cd` command  
followed by the tilde `~`. Try it out, copy paste the following on the  
command line to change to your home directory and list `ls` the files in  
your home environment:  



```r
cd ~  
ls  
```

Your home directory has 1Tb of storage provided free of charge when you open an account  
on MSU HPCC. To check how much space you have available use the disk free  
command `df` with `-h` which prints out the system disk space usage  
in human readble form. Try it out, copy paste the following comand:  


```r
df -h ~
```

  
### Create new directory  
  
Now lets create a new directory. As of now you do not have any directories in your home space.  
You can check this by using the list command `ls`. Try it out:  


```r
ls
```

Think of the name you will give your new directory. Always avoid using spaces in your directories   
and file names. It complicates calling your directories and files. Instead of using a space you  
can use the underscore `_` or hyphen `-`. To create a new directory use the make directory command `mkdir`  
followed by the name of your new directory.  


```r
mkdir name_of_choice 
```

Great job! You just created your first directory in your home space.  
Check out your new directory by using the list command `ls` you already learned.  
By adding `-ltrh` to the `ls` command you are able to see more information on each  
file like permissions, owner, group, size, date modified and filename. Try it out:  


```r
ls -ltrh ~/name_of_choice
```

To go into your newly created folder, similar to a double click on your mouse,  
use the change directory command `cd`  


```r
cd ~/name_of_choice
```

Now you are in your new directory. To see the path to this directory use the print  
working directory command `pwd`. Just write 'pwd' on the command line and it will show  
you the path of the current directory:  


```r
pwd
```

  
### Creating Files  
  
You have now learned how to check the available space on your home directory with `df -h`,  
create a new directory with `mkdir`, check the path to that directory with `pwd` and change  
to a different directory with `cd`. Now we will learn how to create a file in your new directory  
and read it. There are several ways to create a new file from the command line. We will start off  
with one of the most basic ways by using the `nano` command which opens up an editor window.  
To try it, write `nano` on the command line followed by the name of your new file.  
This will open an editor window. Once in, write a few lines of text. To save your work  
press `Ctrl o` and if you have not written your file name when running the `nano` command  
you can write the name of your file. For example: my_text.txt   
Press `enter` to save your text to the file and return to the editor. To exit press `Ctrl x`.  


```r
nano my_text.txt
```

Great! You just saved a file to your new directory in your home space.  
To read your file you can use the `cat` command, which reads your file to standard output.  
Try it by running the following command followed by the name of your file.  


```r
cat my_text.txt
```

There are time when a file contains thousand to million of lines, like in fastq sequence files.  
While their typically is no good reason to try to read a file with thousands of lines on the shell  
you might want to look at a large file from time to time. In these cases you should not use  
the 'cat' command because it will print the entire file to standard output. To avoid this you  
should use the `less` command instead of the `cat` command. Try it out. Copy and paste the following  
code block to your command line. It will modified the file you created.  


```r
fl=(`ls *`)
for ((i=0; i<1000 ; i++ )) do cat $fl >> .tmp; done; mv .tmp $fl 
cat $fl
```

Now that is a lot of test! This file contained only a few thousand lines. Imagine if it had a million lines.  
You would end up seeing a ton of text thrown your way. Now lets try using the `less` command.  
With `less` you are able to see the file in chunks and scroll through the file. You can exit any  
time by pressing the q button. Lets try it:  


```r
less my_txt.txt 
```

**Side note** If you find yourself using `cat` on an annoyingly large file press `Ctrl C` to stop the  
command. This works for any command you wish to terminate on the shell.  
  
  
Now before we finish I want you to learn one of the most dangerous commands to use. The remove command 'rm'.  
The `rm` command deletes directories or files. To delete a file use:  


```r
rm file_name
```

To delete a directory you need to run the remove command recursively by adding the -R.  
For example to delete a directory use:  


```r
rm -R directory_name  
```

Excellent! This concludes our first tutorial. Before you are done let's go to our research space on HPC.  
This research space holds the scripts and output from the analysis Sudeep and myself have done for the group.  
You can look through the directories and files but PLEASE do not delete any file. If you do by mistake delete an  
important file because it does happen from time to time breathe easy. HPC runs a  
full backup every couple of hours so the files can be recovered for the research and home space. The catch is  
you do need to create a ticket with HPC customer service and it can take a couple of days to recover your files.  
I also keep a backup of the research space on our local hard drives, just in case.   
To go to our research directory use:  


```r
cd /mnt/research/NMDL
```

Hope you enjoyed this tutorial. Any comments or suggestions can be set to velezdeb@msu.edu.  
