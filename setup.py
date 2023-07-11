from setuptools import setup, find_namespace_packages


setup(
    name='auto-irida-upload',
    version='0.1.0-alpha',
    packages=find_namespace_packages(),
    entry_points={
        "console_scripts": [
            "auto-irida-upload = auto_irida_upload.__main__:main",
        ]
    },
    scripts=[],
    package_data={
    },
    install_requires=[
    ],
    description='Automated upload of sequence data to IRIDA platform',
    url='https://github.com/BCCDC-PHL/auto-irida-upload',
    author='Dan Fornika',
    author_email='dan.fornika@bccdc.ca',
    include_package_data=True,
    keywords=[],
    zip_safe=False
)
