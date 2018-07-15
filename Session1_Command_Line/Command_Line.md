# Introduction to Command Line
(images and lesson adapted from http://swcarpentry.github.io/)

## 1. Log into guest account 

```
Asia-MacBook-Pro-2:~ arcova$ ssh guestX@rush.ccr.buffalo.edu
Asia-MacBook-Pro-2:~ arcova$ Password: 
```

## 2. Command line prompt

The command ls gives the contents of your current directory
```
Asia-MacBook-Pro-2:~ arcova$ ls 
NEX2_22_S2_R1_fastqc					
NEX2_22_S2_R1_fastqc.zip				
NEX2_22_S2_R2_fastqc					
NEX2_22_S2_R2_fastqc.zip				
NEX2_23_S3_R1_fastqc					
NEX2_23_S3_R1_fastqc.zip
```
What happens at the prompt when you modify the ls command with the flag "-l" ?

```
Asia-MacBook-Pro-2:~ arcova$ ls -l
```

What happens at the prompt when you modify the ls command with the flag "-l -h" ?


```
Asia-MacBook-Pro-2:~ arcova$ ls -l -h 
```

Getting help for commands 
```
Asia-MacBook-Pro-2:~ arcova$ ls --help
Asia-MacBook-Pro-2:~ arcova$ man ls 
```
What are the other available ls flags?

---

## 3. Moving around directories and navigating files 
![alt text](filesystem.svg)

What happens in typing the command pwd?
```
Asia-MacBook-Pro-2:~ arcova$ pwd
Asia-MacBook-Pro-2:~ arcova$ ls
```
Pick a directory name and type cd -chosen directory name-

What happens after you type this ?

```
Asia-MacBook-Pro-2:~ arcova$ cd <pick a directory name>
```
What happens after you type cd ~/? 

```
Asia-MacBook-Pro-2:~ arcova$ cd ~/
```

Hint, the root directory in *nix systems is symbolized by "/," and at the front of a file or directory name it signifies the root directory is signified. 

When "/" appears inside a name, it serves as a separator, eg. "/Users/arcova/data."

The character ~ (tilde) at the start of a path means “the current user’s home directory”. If my home directory is /Users/arcova, then ~/data is equivalent to /Users/arcova/data.

How would you get a file listing of a directory other than your current directory? 

Try moving into one of the subdirectories of your home directory. 

This command will bring you up one directory level. 
```
Asia-MacBook-Pro-2:~ arcova$ cd ..
```

---

## 4. Working with files and directories 

```
Asia-MacBook-Pro-2:~ arcova$ mkdir test_dir
```

This command will make a new directory called test_dir.
Move into your newly made directory. 

```
Asia-MacBook-Pro-2:~ arcova$ nano my_first_text.txt
```
this command will open a new file called "my\_first\_text.txt" in the text editor nano 
write "Hello World!". 
To save the file press Ctrl-O and return and to exit nano press Ctrl-X.
Run the commands to check your current directory and generate a listing of current files. 

On of the functions of the cat command is to display the contents of the file to stdout (i.e. the screen).
```
Asia-MacBook-Pro-2:~ arcova$ cat my_first_text.txt 
```
This command will make a new copy of "my\_first\_text.txt" with a different file name

```
Asia-MacBook-Pro-2:~ arcova$ cp my_first_text.txt my_other_text.txt
```

This command will remove the file "my\_first\_text.txt"

```
Asia-MacBook-Pro-2:~ arcova$ rm my_first_text.txt
Asia-MacBook-Pro-2:~ arcova$  ls
```

## 5. A simple bioinformatics command 

```
Asia-MacBook-Pro-2:~ arcova$ fastqc -o fastqc_results example.fastq.gz
```
a. fastqc indicates which program to run 

b. -o indicates the prefix for output files and fastqc_results is the desired prefix

c. example.fastq.gz indicates the file to process by fastqc

---







