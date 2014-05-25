/*
Copyright 2010-2011, D. E. Shaw Research.
All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are
met:

* Redistributions of source code must retain the above copyright
  notice, this list of conditions, and the following disclaimer.

* Redistributions in binary form must reproduce the above copyright
  notice, this list of conditions, and the following disclaimer in the
  documentation and/or other materials provided with the distribution.

* Neither the name of D. E. Shaw Research nor the names of its
  contributors may be used to endorse or promote products derived from
  this software without specific prior written permission.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE,
DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY
THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT
(INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*/
#ifndef UTIL_OPENCL_H__
#define UTIL_OPENCL_H__
/*
 * has a couple of utility functions to setup and teardown OpenCL.
 * Avoid much boilerplate in every OpenCL program
 */

#if defined(__APPLE__) || defined(__MACOSX)
#include <OpenCL/cl.h>
#else
#include <CL/cl.h>
#endif

#include <cmath>
#include "util.h"
#include "util_macros.h"

#define UCL_STRSIZE 128

// The work group size used by default
#define WG_SIZE 128

typedef struct ucl_info {
    cl_context ctx;
    cl_program prog;
    cl_device_id devid;
    cl_command_queue cmdq;
    cl_uint clkfreq, compunits;
    size_t wgsize;
    cl_ulong localmem;    //local memory size in bytes
    int cores;
    double cycles;
    char vendor[UCL_STRSIZE], devname[UCL_STRSIZE],
	version[UCL_STRSIZE], driver[UCL_STRSIZE];
    int computeflags;
    int rank;
} UCLInfo;

/* Miscellaneous checking macros for convenience */
static char *print_cl_errstring(cl_int err) {
    switch (err) {
        case CL_SUCCESS:                          return strdup("Success!");
        case CL_DEVICE_NOT_FOUND:                 return strdup("Device not found.");
        case CL_DEVICE_NOT_AVAILABLE:             return strdup("Device not available");
        case CL_COMPILER_NOT_AVAILABLE:           return strdup("Compiler not available");
        case CL_MEM_OBJECT_ALLOCATION_FAILURE:    return strdup("Memory object allocation failure");
        case CL_OUT_OF_RESOURCES:                 return strdup("Out of resources");
        case CL_OUT_OF_HOST_MEMORY:               return strdup("Out of host memory");
        case CL_PROFILING_INFO_NOT_AVAILABLE:     return strdup("Profiling information not available");
        case CL_MEM_COPY_OVERLAP:                 return strdup("Memory copy overlap");
        case CL_IMAGE_FORMAT_MISMATCH:            return strdup("Image format mismatch");
        case CL_IMAGE_FORMAT_NOT_SUPPORTED:       return strdup("Image format not supported");
        case CL_BUILD_PROGRAM_FAILURE:            return strdup("Program build failure");
        case CL_MAP_FAILURE:                      return strdup("Map failure");
        case CL_INVALID_VALUE:                    return strdup("Invalid value");
        case CL_INVALID_DEVICE_TYPE:              return strdup("Invalid device type");
        case CL_INVALID_PLATFORM:                 return strdup("Invalid platform");
        case CL_INVALID_DEVICE:                   return strdup("Invalid device");
        case CL_INVALID_CONTEXT:                  return strdup("Invalid context");
        case CL_INVALID_QUEUE_PROPERTIES:         return strdup("Invalid queue properties");
        case CL_INVALID_COMMAND_QUEUE:            return strdup("Invalid command queue");
        case CL_INVALID_HOST_PTR:                 return strdup("Invalid host pointer");
        case CL_INVALID_MEM_OBJECT:               return strdup("Invalid memory object");
        case CL_INVALID_IMAGE_FORMAT_DESCRIPTOR:  return strdup("Invalid image format descriptor");
        case CL_INVALID_IMAGE_SIZE:               return strdup("Invalid image size");
        case CL_INVALID_SAMPLER:                  return strdup("Invalid sampler");
        case CL_INVALID_BINARY:                   return strdup("Invalid binary");
        case CL_INVALID_BUILD_OPTIONS:            return strdup("Invalid build options");
        case CL_INVALID_PROGRAM:                  return strdup("Invalid program");
        case CL_INVALID_PROGRAM_EXECUTABLE:       return strdup("Invalid program executable");
        case CL_INVALID_KERNEL_NAME:              return strdup("Invalid kernel name");
        case CL_INVALID_KERNEL_DEFINITION:        return strdup("Invalid kernel definition");
        case CL_INVALID_KERNEL:                   return strdup("Invalid kernel");
        case CL_INVALID_ARG_INDEX:                return strdup("Invalid argument index");
        case CL_INVALID_ARG_VALUE:                return strdup("Invalid argument value");
        case CL_INVALID_ARG_SIZE:                 return strdup("Invalid argument size");
        case CL_INVALID_KERNEL_ARGS:              return strdup("Invalid kernel arguments");
        case CL_INVALID_WORK_DIMENSION:           return strdup("Invalid work dimension");
        case CL_INVALID_WORK_GROUP_SIZE:          return strdup("Invalid work group size");
        case CL_INVALID_WORK_ITEM_SIZE:           return strdup("Invalid work item size");
        case CL_INVALID_GLOBAL_OFFSET:            return strdup("Invalid global offset");
        case CL_INVALID_EVENT_WAIT_LIST:          return strdup("Invalid event wait list");
        case CL_INVALID_EVENT:                    return strdup("Invalid event");
        case CL_INVALID_OPERATION:                return strdup("Invalid operation");
        case CL_INVALID_GL_OBJECT:                return strdup("Invalid OpenGL object");
        case CL_INVALID_BUFFER_SIZE:              return strdup("Invalid buffer size");
        case CL_INVALID_MIP_LEVEL:                return strdup("Invalid mip-map level");
        default:                                  return strdup("Unknown");
    }
} 

#define CHECKERR(x) do { \
    (x); \
    if (err != CL_SUCCESS) { \
	fprintf(stderr, "Error %d: %s from %s\n", err, print_cl_errstring(err), #x); \
	exit(1); \
    } \
} while(0)

#define CHECK(x) CHECKERR(err = (x))

static UCLInfo *opencl_init(const char *devstr, const char *src, const char *options, cl_uint *ndevices, int debug = 0)
{
#define UCL_MAX_PROPERTIES 32
#define UCL_MAX_PLATFORMS 8
#define UCL_MAX_DEVICES 8
    //UCLInfo *tp;
    UCLInfo *deviceInfo;
    cl_context_properties ctxprop[UCL_MAX_PROPERTIES];
    cl_int err;
    cl_platform_id platforms[UCL_MAX_PLATFORMS];
    cl_uint nplatforms;
    cl_device_id devices[UCL_MAX_DEVICES];
    const char *srcstr[1];
    int cores, devcores;
    double minCycles = 0;

    /* get list of platforms */
    CHECK(clGetPlatformIDs(UCL_MAX_PLATFORMS, platforms, &nplatforms));
    
    if (nplatforms == 0) {
        fprintf(stderr, "No OpenCL platforms available\n");
        return NULL;
    }
    dprintf(("found %d platform%s:\n", nplatforms, nplatforms == 1 ? "" : "s"), debug);
    
    //CHECKNOTZERO(tp = (UCLInfo *)malloc(sizeof(UCLInfo)));

    ctxprop[0] = CL_CONTEXT_PLATFORM;
    ctxprop[1] = 0; /* will fill in platform in loop */
    ctxprop[2] = 0;
    cores = devcores = 0;
    char openCLversion[UCL_STRSIZE];
    int validDevices = 0;
    
    for (int i = 0; i < nplatforms; i++) {
        dprintf(("platform %d: 0x%lx\n", i, (unsigned long)platforms[i]), debug);
        CHECK(clGetDeviceIDs(platforms[i], CL_DEVICE_TYPE_ALL, UCL_MAX_DEVICES, devices, ndevices));
        dprintf(("platform 0x%lx has %d devices:\n", (unsigned long)platforms[i], *ndevices), debug);
        
        ctxprop[1] = (cl_context_properties) platforms[i];
        
        
        for(int j = 0; j < *ndevices; j++) {
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, sizeof openCLversion, openCLversion, 0));
            if(strstr(openCLversion, "OpenCL 1.0") == NULL) {
                validDevices++;
            }
        }
        std::cout << "Start valid devices : " << validDevices << std::endl;
        
        deviceInfo = (UCLInfo *)malloc(sizeof(UCLInfo)*validDevices);
        
        int validDeviceCounter = 0;
        for (int j = 0; j < *ndevices; j++) {
            UCLInfo uc;
            uc.devid = devices[j];
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_NAME, sizeof uc.devname, uc.devname, 0));
            nameclean(uc.devname);
            
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_VENDOR, sizeof uc.vendor, uc.vendor, 0));
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_VERSION, sizeof uc.version, uc.version, 0));
            CHECK(clGetDeviceInfo(devices[j], CL_DRIVER_VERSION, sizeof uc.driver, uc.driver, 0));
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_MAX_CLOCK_FREQUENCY, sizeof uc.clkfreq, &uc.clkfreq, 0));
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_MAX_COMPUTE_UNITS, sizeof uc.compunits, &uc.compunits, 0));
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_MAX_WORK_GROUP_SIZE, sizeof uc.wgsize, &uc.wgsize, 0));
            CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_LOCAL_MEM_SIZE, sizeof uc.localmem, &uc.localmem, 0));
            uc.computeflags = 0;
            uc.wgsize = WG_SIZE;
            
#ifdef CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV
            {
                cl_uint nvmaj, nvmin;
                CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_COMPUTE_CAPABILITY_MAJOR_NV, sizeof nvmaj, &nvmaj, 0));
                CHECK(clGetDeviceInfo(devices[j], CL_DEVICE_COMPUTE_CAPABILITY_MINOR_NV, sizeof nvmin, &nvmin, 0));
                uc.computeflags = nvmaj*10 + nvmin;
            }
#endif
            cores = uc.compunits;
            /* XXX Hardwired knowledge about devices */
            if (strcmp("Cayman", uc.devname) == 0) {
                /*
                 * Each modern AMD compute unit (shader cluster?)
                 * is a 16-lane SIMD Engine with 4 (Cayman) or 5
                 * (Cypress) VLIW slots per lane.  AMD appears to
                 * think of each slot as a "stream processor"
                 * (shader processor) in their marketing i.e. a
                 * Cayman-based Radeon 6950 with 24 compute units
                 * has 1536 stream processors.
                 */
                cores *= 16*4;
            } else if (strcmp("Cypress", uc.devname) == 0) {
                cores *= 16*5;
            } else if (strstr(uc.devname, "GTX 480") ||
                   strstr(uc.devname, "C20") ||
                   strstr(uc.devname, "M20")) {
                /*
                 * Fermi has 32 cores per SM.  Maybe use
                 * computeflags to figure this out?
                 */
                cores *= 32;
            } else if(strstr(uc.devname, "GTX 580")) {
                // Manually do fixes for GTX 580 on Mac Lion
                uc.compunits = 16;
                cores = 16*32;
                uc.clkfreq = 1564;       // ASUS GTX 580 is overclocked to 782 MHz by default this is showing up as 0 MHz
                
            } else if(strstr(uc.devname, "GTX 680")) {
                uc.compunits = 8;
                cores = 8*192;  // = streaming multiprocessors (SM's x Number of Cores per SM)
                uc.clkfreq = 1058;
            } else if(strstr(uc.devname, "Intel")) {
                uc.wgsize = 1;
            }
            
            
            /* clkfreq is in Megahertz! */
            uc.cycles = 1e6 * uc.clkfreq * cores;
            dprintf(("  %d: device 0x%lx vendor %s %s version %s driver %s : %u compute units @ %u MHz %d cores cycles/s %.2f flags %d\n", j, (unsigned long) devices[j], uc.vendor, uc.devname, uc.version, uc.driver, uc.compunits, uc.clkfreq, cores, uc.cycles, uc.computeflags), debug);
            /*
            printf("Device ID : 0x%lx : %s : %d units %d cores %.2f Gcycles/s, %lu max work group size %lu bytes local memory\n", (unsigned long)tp->devid, tp->devname, tp->compunits, devcores, tp->cycles*1e-9, tp->wgsize, tp->localmem);
            */
            if (devstr && strstr(uc.devname, devstr) == NULL) {
                if (debug) {
                    printf("skipping device %s\n", uc.devname);
                }
                continue;
            }
            
            if(strstr(uc.version, "OpenCL 1.0") != NULL) {
                std::cout << " Found a bad one " << std::endl;
                continue;
            } else {
            
                deviceInfo[validDeviceCounter].devid = uc.devid;
                deviceInfo[validDeviceCounter].clkfreq = uc.clkfreq;
                deviceInfo[validDeviceCounter].compunits = uc.compunits;
                deviceInfo[validDeviceCounter].wgsize = uc.wgsize;
                deviceInfo[validDeviceCounter].localmem = uc.localmem;
                deviceInfo[validDeviceCounter].cores = uc.cores;
                deviceInfo[validDeviceCounter].cycles = uc.cycles;
                deviceInfo[validDeviceCounter].localmem = uc.localmem;
                strcpy(deviceInfo[validDeviceCounter].vendor, uc.vendor);
                strcpy(deviceInfo[validDeviceCounter].devname, uc.devname);
                strcpy(deviceInfo[validDeviceCounter].version, uc.version);
                strcpy(deviceInfo[validDeviceCounter].driver, uc.driver);
                deviceInfo[validDeviceCounter].computeflags = uc.computeflags;
            
                // Create a context
                CHECKERR(deviceInfo[j].ctx = clCreateContext(ctxprop, 1, &uc.devid, 0, 0, &err));
                dprintf(("create OpenCL context for device 0x%lx %s\n", (unsigned long)deviceInfo[j].devid, deviceInfo[j].devname), debug);
                
                // Create a command queue
                CHECKERR(deviceInfo[j].cmdq = clCreateCommandQueue(deviceInfo[j].ctx, deviceInfo[j].devid, 0, &err));
                
                // Build program for device
                dprintf(("create OpenCL program from source\n"), debug);
                srcstr[0] = src;
                CHECKERR(deviceInfo[j].prog = clCreateProgramWithSource(deviceInfo[j].ctx, 1, srcstr, 0, &err));
                fprintf(stderr, "Create Program Err : %d\n",err);
                //CHECK(clBuildProgram(tp->prog, 1, &tp->devid, options, 0, 0));
                
                if (clBuildProgram(deviceInfo[j].prog, 1, &deviceInfo[j].devid, options, 0, 0) != CL_SUCCESS) {
                    char errbuf[512*1024];
                    cl_int builderr = err;
                    CHECK(clGetProgramBuildInfo(deviceInfo[j].prog, deviceInfo[j].devid, CL_PROGRAM_BUILD_LOG, sizeof errbuf, &errbuf, 0));
                    
                    if (errbuf[0]) {
                        fprintf(stderr, "OpenCL build for device id 0x%lx %s failed with error %d: %s\n", (unsigned long) deviceInfo[j].devid, deviceInfo[j].devname, builderr, errbuf);
                    }
                    
                    CHECK(clGetProgramBuildInfo(deviceInfo[j].prog, deviceInfo[j].devid, CL_PROGRAM_BUILD_LOG, sizeof errbuf, &errbuf, 0));
                    
                    
                    exit(1);
                }
                
                
                // Record the lowest power device
                if(minCycles == 0 || deviceInfo[j].cycles < minCycles) {
                    minCycles = deviceInfo[j].cycles;
                }
            
                validDeviceCounter++;
                std::cout << " Device counter : " << validDeviceCounter << std::endl;
            }
        } // for ndevices
    } // For nplatforms
        
    // Rank the devices
    std::cout << "Number of valid Devices : " << validDevices << std::endl;
    for(int i = 0; i < validDevices; i++) {
        float x = deviceInfo[i].cycles/(float)minCycles;
        
        deviceInfo[i].rank = (x > (floor(x)+0.5f)) ? ceil(x) : floor(x);
        printf("Device %s rank %d\n", deviceInfo[i].devname, deviceInfo[i].rank);
    }
    
    *ndevices = validDevices;
    
    dprintf(("opencl_init done\n"), debug);
    /* XXX Save build programs as .deviceid so we can read them back and run? */
    return deviceInfo;
}


static void opencl_done(UCLInfo *tp, cl_uint ndevices, int debug = 0) {
    cl_int err;
    
    dprintf(("opencl_done\n"), debug);
    
    for(int  i = 0; i < ndevices; i++) {
        CHECK(clReleaseCommandQueue(tp[i].cmdq));
        tp[i].cmdq = 0;
        CHECK(clReleaseProgram(tp[i].prog));
        tp[i].prog = 0;
        CHECK(clReleaseContext(tp[i].ctx));
        tp[i].ctx = 0;
    }
    
    free(tp);
}


#endif /* UTIL_OPENCL_H__ */
