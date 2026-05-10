from setuptools import setup, find_packages

setup(
    name="AlperenMSA", 
    version="0.1.0",
    author="Adın Soyadın",
    author_email="mailadresin@email.com",
    description="Dinamik Programlama ile Çoklu Dizi Hizalaması (MSA) paketi",
    long_description="Bu kütüphane biyoinformatik bahar dönemi projesi kapsamında Needleman-Wunsch ve 3-Dizi DP hizalaması yapar.",
    long_description_content_type="text/markdown",
    packages=find_packages(), 
    install_requires=[
        "numpy"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.6',
)