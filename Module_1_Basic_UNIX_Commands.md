---
title: Basic UNIX Commands
author: Deborah Velez-Irizarry
date: Updated Jan 2 2020
output:
  prettydoc::html_pretty:
    theme: hpstr
    highlight: github
    toc: true
---


### Description

This tutorial is a review of some basic LINUX/UNIX commands, such as creating directories, saving and reading files from the command line, and some useful commands to keep in mind. These are core skills you will use every time you work on HPC.


### Connecting to the HPC system

There are various ways to connect to HPC.

> Web-based remote desktop environment: [Web site access to HPCC](https://wiki.hpcc.msu.edu/display/ITH/Web+Site+Access+to+HPCC)

> SSH client: [iCER](https://wiki.hpcc.msu.edu/display/ITH/Connect+to+HPC+System)

Take some time to learn more about MSU HPCC. Some good resources to look over are:

> [Why use HPCC](https://wiki.hpcc.msu.edu/display/ITH/Why+Use+HPCC)
> [HPC Glossary](https://wiki.hpcc.msu.edu/display/ITH/HPC+Glossary)
> [HPC entire layout at iCER](https://wiki.hpcc.msu.edu/display/ITH/HPC%27s+entire+layout+at+iCER)
> [Training Resources](https://icer.msu.edu/education-events/training-resources)


Since we will be working on the terminal, let us use the SSH client to log in. Once you got Ubuntu up and running, press `Ctrl`+`Alt`+`T` to open the terminal. From the command line, log in to HPCC using `ssh`. They will ask for your MSU password.


```bash
ssh -YX username@gateway.hpcc.msu.edu
```

**Let us get started.**


### Home Directory

Let us start by looking at your home directory on HPC. Logging in to HPC takes you directly to your home directory. The tilde `~` symbol is an alias for your home environment. You can go to your home directory by just using the change directory `cd` command, followed by the tilde `~`. Try it out, copy-paste the following on the command line to change to your home directory and list `ls` the files in your home environment:


```bash
cd ~
ls
```

Your home directory has 1Tb of storage provided free of charge when you open an account on MSU HPCC. To check how much space you have available, use the disk free command `df` with `-h`, which prints out the system disk space usage in a human-readable form. Try it out, copy-paste the following command:


```bash
df -h ~
```


### Create new directory


Think of the name you will give your new directory. Always avoid using spaces in your folder and file names. It complicates calling your directories and files. Instead of using spaces, you can use the underscore `_` or hyphen `-`. To create a new directory, use the make directory command `mkdir` followed by the name of your new directory.


```bash
mkdir name_of_choice
```

Great job! You just created your first directory in your home space. Check out your new directory by using the list command `ls` you already learned. By adding `-ltrh` to the `ls` command, you can see more information on each file like permissions, owner, group, size, date modified, and filename. Try it out:


```bash
ls -ltrh ~/name_of_choice
```

To go into your newly created folder, similar to a double click on your mouse, use the change directory command `cd`.


```bash
cd ~/name_of_choice
```

Now you are in your new directory. To see the path to this directory, use the print working directory command `pwd`. Just write 'pwd' on the command line, and it will show you the path of the current directory:


```bash
pwd
```


### Creating Files

You have now learned how to check the available space on your home directory with `df -h`, create a new directory with `mkdir`, check the path to that directory with `pwd` and change to a different directory with `cd`. Now we will learn how to create a file in your new directory and read it. There are several ways to create a new file from the command line. We will start with one of the most basic means by using the `nano` command, which opens up an editor window. To try it, write `nano` on the command line, followed by the name of your new file. The `nano` command will open an editor window. Once in, write a few lines of text. To save your work press `Ctrl o` and if you have not written your file name when running the `nano` command, you can write the name of your file. For example my_text.txt

```bash
nano my_text.txt
```

> Press `enter` to save your text to the file and return to the editor. To exit press `Ctrl x`.

Great! You just saved a file to your new directory in your home space. To open a file, you can use the `cat` command, which reads your file to standard output. Let us try it, write `cat` on the command line followed by your file name with the extension.

```bash
cat my_text.txt
```

There are times when a file contains thousand to millions of lines, like in fastq sequence files. While there typically is no good reason to try to read a large file directly on the shell, you might want to preview one from time to time. In these cases, you should not use the 'cat' command because it will print the entire file to standard output. To avoid this, you should use the `less` command instead of the `cat` command. Try it out. Copy and paste the following code block to your command line. It will modify the file you created.


```bash
fl=(`ls *`)
for ((i=0; i<1000 ; i++ )) do cat $fl >> .tmp; done; mv .tmp $fl
cat $fl
```

Now that is a lot of text! This file contained only a few thousand lines. Imagine if it had a million lines. You would end up seeing a ton of text thrown your way. Now let us try using the `less` command. With `less` you can see the file in chunks and scroll through the file. You can exit any time by pressing the q button. Let us try it:


```bash
less my_txt.txt
```

**Side note** If you find yourself using `cat` on an annoyingly large file, press `Ctrl C` to stop the command. Use `Ctrl C` on any command you wish to terminate on the shell.


Now before we finish, I want you to learn one of the most dangerous commands to use. The remove command 'rm'. The `rm` command deletes directories or files. To delete a file use:


```bash
rm file_name
```

To delete a directory, you need to run the remove command recursively by adding the -R. For example to eliminate a directory use:


```bash
rm -R directory_name
```

Excellent! Before we finish, let us go to our research space on HPC. This research space holds the scripts and output from the analysis Sudeep and I have done for the group. You can look through the directories and files, but PLEASE do not delete any file. If you happen to delete a file by mistake, it does happen from time to time, breathe easy HPC runs a full backup every couple of hours, so the files are recoverable for the research and home space. The catch is you do need to create a ticket with HPC customer service, and it can take a couple of days to recover your files. I also keep a backup of the research space on our local hard drives, just in case. To go to our research directory use:


```bash
cd /mnt/research/NMDL
```

I hope you enjoyed this tutorial. Send any comments or suggestions to velezdeb@msu.edu.



