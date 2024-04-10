# Installation
Installation should be straightforward. The installation was never tested on Mac but it should work as well.
#### 1 - Install python 
On Windows, we recommend not to use Anaconda but rather to [download python](https://www.python.org/downloads/windows/) and install it manually. Get Python3.9 or newer.
#### 2- Create a new environment
*Why should you create a python environment ?*

![](https://imgs.xkcd.com/comics/python_environment.png)

*Image from [xkcd.com](https://xkcd.com/1987/) by Randall Munroe.*

Because you do not want to be able to break your environment without breaking your machine. Or not to know *where* you installed all your stuffs. So, yes, please, do it.


Open a terminal (Linux) or a shell (Windows) and enter
```
python3 -m  pip install virtualenv
```
Now you can create a new environment. In the following, we will call it *helenv*.
```
python -m venv helenv
```
Will create in your current directory a folder named 'helenv' in which the source code of your environment will be. Now, you need to activate your environment; this will be necessary every time you want to use your Python because we will install all the packages inside this environment.
```
Linux : $ source helenv/bin/activate
Windows : $ helenv\bin\activate
```
Once your environment is activated, you should see it in between parenthesis in you shell.

#### 3- Download heliumtools
The source code can be downloaded from the [gitlab](https://gitlab.in2p3.fr/gaz-quantiques-lcf/helium-1/heliumtools) or the [github](https://github.com/quantumatomoptic/heliumtools) repository. If you do not want to participate to the code, you are not forced to register your SSH key and/or to have a github-gitlab account. Otherwise, you  should register a [SSH key](#install-an-ssh-key).  
#### 4- Install heliumtools and PyTorch
Once you have the code, install it as a developer. This will create symbolic links from your environment to the heliumtools folder and it means that if you pull the repository, it will automatically update your repository. 
```
cd heliumtools 
pip install -e .
pip install torch --index-url https://download.pytorch.org/whl/cpu
```
The `-e` option stands for developer mode and the dot `.` means *here*. Last but not least, please also install PyTorch. The link here assume you have a normal CPU and not a CUDA based one.



## How to section
#### Install an SSH Key
Open a terminal or a cmd window and enter
```
ssh-keygen
```
The default options will create SSH-RSA files a .ssh folder in your HOME (`C:\Users\fname\.ssh` for Windows and `~/.ssh` for Linux). You can then go in tht folder and collect the public key in the file `id_rsa.pub` (open in with a text app). You will then copy the content of this file. Go to gitlab/github, and in your account parameter, add an ssh key. Register the public key you just copied and save it. You can now pull your code to the repository !









