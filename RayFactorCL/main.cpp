//
//  main.cpp
//  RayFactor OpenCl
//
//  Created by Trevor Walker on 1/03/12.
//  Copyright 2012 Native Dynamics. All rights reserved.
//

#include <iostream>
#include <assert.h>
#include <time.h>       // For the gettime
#include <sys/stat.h>   // For the load program source
#include <fstream>      // for printing to file

#if defined(__APPLE__) || defined(__MACOSX)
    #include <OpenCL/cl.h>
#else
    #include <CL/cl.h>
#endif

#include "cycle.h"
#include "util_opencl.h"
#include "Scene.h"

int runCL_GPURNG(int size, Scene &scene);
char * load_program_source(const char *filename);


#pragma mark -
#pragma mark Utilities
char * load_program_source(const char *filename)
{
    struct stat statbuf;
	FILE *fh;
	char *source;
	
	fh = fopen(filename, "r");
	if (fh == 0)
		return 0;
	
	stat(filename, &statbuf);
	source = (char *) malloc(statbuf.st_size + 1);
	fread(source, statbuf.st_size, 1, fh);
	source[statbuf.st_size] = '\0';
	
	return source;
}



#pragma mark -
#pragma mark Main OpenCL Routine

int runCL_GPURNG(int n_samples, Scene &scene)
{
    
    Primitive *objects = &scene.primitives[0];       // Note the requirement that vectors be continuous in memory was included in C++03 (C++98-TC1) 23.2.4.1
    BoundingVolume *bvs = &scene.boundingVolumes[0];
    int n_objects = (int)scene.primitives.size();
    int n_bvs = (int)scene.boundingVolumes.size();
    
    UCLInfo *devices;
    cl_uint ndevices = 0;
    const char * options = "-I ./ -cl-no-signed-zeros -cl-mad-enable -cl-single-precision-constant -cl-strict-aliasing -cl-unsafe-math-optimizations -cl-finite-math-only";
    const char *kernelname;
    const char *filename = "RayFactorGPU.cl";
    
    if(scene.useBoundingVolumes)
    {
        kernelname = "runBVGPU";
    }
    else
    {
        kernelname = "runGPU";
    }
    
    char *program_source = load_program_source(filename);
    
    cl_int err = 0;
    size_t nthreads, buffer_size;
    cl_mem objects_mem, bvs_mem;
    
    ticks t0, t1;
    
    double *samplesForObj = new double[n_objects];
    
    std::cout << "=========== RayFactor CL =============" << std::endl;
    t0 = getticks();
    devices = opencl_init(NULL, program_source, options, &ndevices, 1);
    
    cl_kernel *kernel = new cl_kernel[ndevices];
    cl_mem *hits_mem = new cl_mem[ndevices];
    
    free(program_source);
    t1 = getticks();
    std::cout << "OpenCL initalise : " << elapsed(t1,t0) << " cycles" << std::endl;
    
    
    // This method you specify the number of iterations each work item does
    
    int count = 4; // Number of elements in the random number
    
    int *hits = new int[ndevices*n_objects*n_objects];
    
    for (int i =0; i < ndevices*n_objects*n_objects; i++)
    {
        hits[i] = 0;
    }
    
    double totalRank = 0;

    // Set the kernel arguments for each device
    for (int i = 0; i < ndevices; i++)
    {
        totalRank += devices[i].rank;
        CHECKERR(kernel[i] = clCreateKernel(devices[i].prog, kernelname, &err));
        
        // Crete the memory buffers
        buffer_size = sizeof(Primitive)*n_objects;
        objects_mem = clCreateBuffer(devices[i].ctx, CL_MEM_READ_ONLY, buffer_size, NULL, NULL);
        err = clEnqueueWriteBuffer(devices[i].cmdq, objects_mem, CL_TRUE, 0, buffer_size, (void*)objects, 0, NULL, NULL);
        
        if(n_bvs > 0)
        {
            buffer_size = sizeof(BoundingVolume)*n_bvs;
            bvs_mem = clCreateBuffer(devices[i].ctx, CL_MEM_READ_ONLY, buffer_size, NULL, NULL);
            err = clEnqueueWriteBuffer(devices[i].cmdq, bvs_mem, CL_TRUE, 0, buffer_size, (void*)bvs, 0, NULL, NULL);
        }
        
        buffer_size = sizeof(int) * n_objects*n_objects;
        CHECKERR(hits_mem[i] = clCreateBuffer(devices[i].ctx, CL_MEM_WRITE_ONLY, buffer_size, NULL, &err));
        
        void *hitsLoc = &(hits[i*n_objects*n_objects]);
        CHECKERR(clEnqueueWriteBuffer(devices[i].cmdq, hits_mem[i], CL_TRUE, 0, buffer_size, (void*)hitsLoc, 0, NULL, NULL));
        
        CHECK(clFinish(devices[i].cmdq));
        
        // Set the args on the device
        CHECK(clSetKernelArg(kernel[i], 1, sizeof(int), &n_objects));
        CHECK(clSetKernelArg(kernel[i], 2, sizeof(cl_mem), &objects_mem));
        CHECK(clSetKernelArg(kernel[i], 3, sizeof(int)*n_objects, NULL));
        CHECK(clSetKernelArg(kernel[i], 4, sizeof(cl_mem), (void *)&(hits_mem[i])));
        
        if(n_bvs > 0)
        {
            CHECK(clSetKernelArg(kernel[i], 5, sizeof(int), &n_bvs));
            CHECK(clSetKernelArg(kernel[i], 6, sizeof(cl_mem), &bvs_mem));
        }
    }
    
    
    // Now setup the arguments to our kernel
    
    std::cout << "Number of objects : " << n_objects << std::endl;
    std::cout << "Number of devices : " << ndevices << std::endl;
    t1 = getticks();
    std::cout << "Load Kernel and memory : " << elapsed(t1,t0) << " cycles" << std::endl;
    
    
    t0 = getticks();
    
    int objNo = 0;
    for(int i = 0; i < ndevices; i++)
    {
        int noObjToProcess = floor(devices[i].rank * n_objects/totalRank + 0.5f);
        
        std::cout << " == Device : " << i << ", Objects to process : " << noObjToProcess<< std::endl;
        for(int j = 0; j < noObjToProcess && objNo < n_objects; j++)
        {
            CHECK(clSetKernelArg(kernel[i], 0, sizeof(int), &objNo));
            
            nthreads = objects[i].rayDensity*n_samples/count;
            nthreads = nthreads % devices[i].wgsize > devices[i].wgsize/2 ? nthreads + devices[i].wgsize - nthreads % devices[i].wgsize : nthreads - nthreads % devices[i].wgsize;
            
            
            CHECK(clEnqueueNDRangeKernel(devices[i].cmdq, kernel[i], 1, 0, &nthreads, &devices[i].wgsize, 0, 0, 0));
            
            samplesForObj[objNo] = nthreads*count;
            
            objNo++;
        }
    }
    
    // Issue finish commands
    for(int i = 0; i < ndevices; i++)
    {
        clFinish(devices[i].cmdq);
    }
    
    t1 = getticks();
    float t = elapsed(t1,t0);
    std::cout << "Enqueue and execute : " << t << std::endl;
    
    
    // Once finished read back sthe results from the answer array into the results array
    
    t0 = getticks();
    
    for (int i = 0; i < ndevices; i++)
    {
        buffer_size = sizeof(int) * n_objects*n_objects;
        CHECK(clEnqueueReadBuffer(devices[i].cmdq, hits_mem[i], CL_TRUE, 0, sizeof(int) * n_objects*n_objects, &hits[i*n_objects*n_objects], 0, NULL, NULL));
        assert(err == CL_SUCCESS);
        clFinish(devices[i].cmdq);
    }
    
    t1 = getticks();
    
    
    std::cout << "Read Memory : " << elapsed(t1,t0) << " cycles" << std::endl;
    
    // Clean up context and free memory
    opencl_done(devices, ndevices);
    
    
    float *vfm = new float[n_objects*n_objects];
    
    for (int i = 0; i < n_objects*n_objects; i++)
    {
        vfm[i] = 0.0f;
    }
    
    // Combine the results from each of the devices
    for (int i = 0; i < ndevices; i++)
    {
        for (int j = 0; j < n_objects; j++)
        {
            for (int w = 0; w < n_objects; w++)
            {
                int index = i*n_objects*n_objects+j*n_objects+w;
                vfm[j*n_objects+w] += hits[index];
            }
        }
    }
    
    std::ofstream out;
    out.open("output.txt");
	out.precision(8);
    
    // Divide by the number of samples used per object and output the result
    for (int i = 0; i < n_objects; i++)
    {
        for (int j = 0; j < n_objects; j++)
        {
            vfm[i*n_objects+j] = vfm[i*n_objects+j]/samplesForObj[i];
            out << vfm[i*n_objects + j] << ",";
        }
        out << std::endl;
    }
    
    out.close();
    
    // Clean up memory
    clReleaseMemObject(objects_mem);
    
    if(n_bvs > 0)
    {
        clReleaseMemObject(bvs_mem);
    }
    
    for (int i = 0; i < ndevices; i++)
    {
        clReleaseMemObject(hits_mem[i]);
    }

    delete hits;
    delete samplesForObj;
    delete kernel;
    delete hits_mem;
    
    return CL_SUCCESS;
}


int main (int argc, const char * argv[])
{
    Scene scene;
	const char *xmlpath = "input.xml";
	scene.readScene(xmlpath);
    
    runCL_GPURNG(scene.globalRayDensity, scene);

    return 0;
}

