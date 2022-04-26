rule compile_mosaic:
    """
    rule fct:
    input:
    output:
    """
    input: 
        "scripts/Mosaic/src/CMakeLists.txt"
    output:
        "scripts/Mosaic/build/mosaic"
    conda:
        "../envs/cpp.yaml"
        # "mc-base"
    shell:
        """
        pwd
        cd scripts/Mosaic/build
        cmake ../src
        make
        cd ../../..
        pwd
        """
