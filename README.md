# IMCprocess
Python package for downstream analysis of Imaging Mass Cytometry files

## IMCprocess docker container
### Requirement 
* Docker installed on your local machine
### Usage
1. Clone or download the IMCprocess repository to your local machine.
2. Build the Docker image using the following command:
```bash
docker build -t imcprocess-jupyter .
```
3. Once the image is built, start a new Docker container using the following command:
```bash
docker run -p 8888:8888 -v /path/to/analysis_dir:/app/analysis_dir -it imcprocess-jupyter
```
Replace /path/to/analysis_dir with the path to the directory on your local machine where you want to store your Jupyter notebooks.
This will start a new Docker container and mount the /path/to/analysis_dir directory on your local machine to the /app/analysis_dir directory inside the container. The -p option maps port 8888 inside the container to port 8888 on your local machine, allowing you to access the Jupyter Notebook server from your web browser.

4. In a web browser on your local machine, navigate to the following URL to access the Jupyter Notebook interface:
`http://localhost:8888/`
5. To create a new Python 3 notebook, click on "New" in the upper right corner and select "Python 3". You can then add the /app directory to your Python path in the Jupyter notebook to import the IMCprocess package:
```python
import sys
sys.path.append('/app')
```

### Note
* All notebooks created inside the container will be saved to the directory specified in the docker run command. To make sure your notebooks are saved on your local machine, be sure to specify a path that is accessible from within the container.
* If you need to install additional packages, you can do so using the pip command inside the container. Any packages installed inside the container will not persist after the container is stopped or removed, unless you save the container state as a new image using the docker commit command.
