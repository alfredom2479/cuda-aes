#include <stdio.h>
#include <stdlib.h>
#include <fcntl.h>
#include <cuda_runtime.h>
#include <time.h>
#include "string.h"

#define DEFAULT_FILENAME "gary.ppm"

#define BLOCK_SIZE 256

unsigned char *read_ppm(char *filename, int *xsize, int *ysize, int *maxval){
    FILE *fp = fopen(filename, "rb");
    if(!fp){
        fprintf(stderr, "Error: '%s' cannot be opened\n", filename);
        return NULL;
    }
    char first_line[4];
    fgets(first_line,4,fp);
    if(strcmp(first_line, "P6\n")){
        fprintf(stderr, "Error: '%s' not in P6 format\n", filename);
        return NULL;
    }
    char temp;
    while((temp = fgetc(fp)) == '#'){
        fscanf(fp, "%*[^\n]\n");
    }
    ungetc(temp,fp);
    fscanf(fp, "%d %d\n%d\n", xsize, ysize, maxval);
    
    unsigned int *pic  = (unsigned int*)malloc((*xsize)*(*ysize)*sizeof(unsigned int));
    unsigned char *ppic = (unsigned char*)malloc((*xsize)*(*ysize)*sizeof(unsigned char)*3);
    for(int i = 0; i < (*xsize)*(*ysize); i++){
        unsigned char buf[3];
        fread((void*)&buf[0], 3, 1, fp);
        //printf("%x,%x,%x\n",buf[0],buf[1],buf[2]);
        //pic[i] = buf[0];
        ppic[i*3]= buf[0];
        ppic[i*3+1] = buf[1];
        ppic[i*3+2] = buf[2];
        
    }
    return ppic;
}

void write_ppm(unsigned char *filename, int xsize, int ysize, int maxval, unsigned char *pic){
    FILE *fp;
    
    fp = fopen((char *)filename, "wb");
    if(!fp){
        fprintf(stderr, "FAILED TO OPEN FILE '%s' for writing\n", filename);
        exit(-1);
    }
    
    fprintf(fp, "P6\n");
    fprintf(fp, "%d %d\n%d\n", xsize, ysize, maxval);
    
    int numpix = xsize * ysize;
    for(int i=0; i < numpix; i++){
        //unsigned char uc = (unsigned char) pic[i];
        fprintf(fp, "%c%c%c", pic[i*3], pic[i*3+1], pic[i*3+2]);
    }
    
    fclose(fp);
}


//__global__ void aes(unsigned int* pic, int* r_pic, int xsize, int ysize, int thresh){
__global__ void aes_naive(unsigned char* pic, unsigned char* r_pic, unsigned char* keys, unsigned char* sbox){
    //current plan is to prepend a 0x00 per pixel so that we can have split up data evenly
    //4 pixels per thread
    //first 16 bytes will contain the number of extra bytes that have been appended to the array
    unsigned char res[4];
    unsigned char mcMatrix[] ={
        0x02, 0x03, 0x01, 0x01,
        0x01, 0x02, 0x03, 0x01,
        0x01, 0x01, 0x02, 0x03,
        0x03, 0x01, 0x01, 0x02
    };
    
    int tx = (blockIdx.x*BLOCK_SIZE+threadIdx.x)*16;
    unsigned char temp, temp0, temp1;
    
    //make local copy to
    unsigned char l_pic[16];
    
    for(int i = 0; i< 16; i++){
        l_pic[i] = pic[tx+i];    
    }
    //initial add round key
    for(int i = 0; i < 16; i++){
        l_pic[i] = l_pic[i] ^ keys[i];
    }
    //The main aes loop
    for(int i = 0; i < 10; i++){
        //sub bytes
        for(int j = 0; j < 16; j++){
            l_pic[j] = sbox[l_pic[j]];
        }
        //shift rows
        temp = l_pic[1];
        l_pic[1] = l_pic[1+4];
        l_pic[1+4] = l_pic[1+8];
        l_pic[1+8] = l_pic[1+12];
        l_pic[1+12] = temp;
        
        temp = l_pic[2];
        l_pic[2] = l_pic[2+8];
        l_pic[2+8] = temp;
        temp = l_pic[2+4];
        l_pic[2+4] = l_pic[2+12];
        l_pic[2+12] = temp;
        
        temp = l_pic[3];
        l_pic[3] = l_pic[3+12];
        l_pic[3+12] = l_pic[3+8];
        l_pic[3+8] = l_pic[3+4];
        l_pic[3+4] = temp;
        
        //mix columns
        if(i<9){
            for(int j=0; j < 4; j++){
                res[0] = 0x00;
                res[1] = 0x00;
                res[2] = 0x00;
                res[3] = 0x00;
                
                for(int k = 0; k < 4; k++){
                    //res[0] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[0+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1>>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[0] = res[0] ^ temp;
                    
                    //res[1] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[4+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[1] = res[1] ^ temp;
                    
                    //res[2] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[8+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[2] = res[2] ^ temp;
                    
                    //res[3] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[12+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[3] = res[3] ^ temp;
                }
                l_pic[(j*4)+0] = res[0];
                l_pic[(j*4)+1] = res[1];
                l_pic[(j*4)+2] = res[2];
                l_pic[(j*4)+3] = res[3];
            }
            
        }
        //add round key
        for(int j = 0; j < 16; j++){
            l_pic[j] = l_pic[j] ^ keys[((i+1)*16)+j];
        }
       //end of of current round of encryption
    }
    //end of all encryption
    //move local bytes back to global memory
    for(int i = 0; i < 16; i++){
        r_pic[tx+i] =  l_pic[i];
    }
}

__global__ void aes(unsigned char* pic, unsigned char* r_pic, unsigned char* keys, unsigned char* sbox){
    //current plan is to prepend a 0x00 per pixel so that we can have split up data evenly
    //4 pixels per thread
    //first 16 bytes will contain the number of extra bytes that have been appended to the array
    
    __shared__ unsigned char keys_s[176];
    __shared__ unsigned char sbox_s[256];
    
    sbox_s[threadIdx.x] = sbox_s[threadIdx.x];
    
    if(threadIdx.x < 176){
        keys_s[threadIdx.x] = keys_s[threadIdx.x];
    }
    unsigned char res[4];
    unsigned char mcMatrix[] ={
        0x02, 0x03, 0x01, 0x01,
        0x01, 0x02, 0x03, 0x01,
        0x01, 0x01, 0x02, 0x03,
        0x03, 0x01, 0x01, 0x02
    };
    
    int tx = (blockIdx.x*BLOCK_SIZE+threadIdx.x)*16;
    unsigned char temp, temp0, temp1;
    
    //make local copy to
    unsigned char l_pic[16];
    
    for(int i = 0; i< 16; i++){
        l_pic[i] = pic[tx+i];    
    }
    //initial add round key
    for(int i = 0; i < 16; i++){
        l_pic[i] = l_pic[i] ^ keys[i];
    }
    //The main aes loop
    for(int i = 0; i < 10; i++){
        //sub bytes
        for(int j = 0; j < 16; j++){
            l_pic[j] = sbox_s[l_pic[j]];
        }
        //shift rows
        temp = l_pic[1];
        l_pic[1] = l_pic[1+4];
        l_pic[1+4] = l_pic[1+8];
        l_pic[1+8] = l_pic[1+12];
        l_pic[1+12] = temp;
        
        temp = l_pic[2];
        l_pic[2] = l_pic[2+8];
        l_pic[2+8] = temp;
        temp = l_pic[2+4];
        l_pic[2+4] = l_pic[2+12];
        l_pic[2+12] = temp;
        
        temp = l_pic[3];
        l_pic[3] = l_pic[3+12];
        l_pic[3+12] = l_pic[3+8];
        l_pic[3+8] = l_pic[3+4];
        l_pic[3+4] = temp;
        
        //mix columns
        if(i<9){
            for(int j=0; j < 4; j++){
                res[0] = 0x00;
                res[1] = 0x00;
                res[2] = 0x00;
                res[3] = 0x00;
                
                for(int k = 0; k < 4; k++){
                    //res[0] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[0+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1>>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[0] = res[0] ^ temp;
                    
                    //res[1] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[4+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[1] = res[1] ^ temp;
                    
                    //res[2] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[8+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[2] = res[2] ^ temp;
                    
                    //res[3] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[12+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[3] = res[3] ^ temp;
                }
                l_pic[(j*4)+0] = res[0];
                l_pic[(j*4)+1] = res[1];
                l_pic[(j*4)+2] = res[2];
                l_pic[(j*4)+3] = res[3];
            }
            
        }
        //add round key
        for(int j = 0; j < 16; j++){
            l_pic[j] = l_pic[j] ^ keys_s[((i+1)*16)+j];
        }
       //end of of current round of encryption
    }
    //end of all encryption
    //move local bytes back to global memory
    for(int i = 0; i < 16; i++){
        r_pic[tx+i] =  l_pic[i];
    }
}


__global__ void aes_moreWork(unsigned char* pic, unsigned char* r_pic, unsigned char* keys, unsigned char* sbox, int size){
    //current plan is to prepend a 0x00 per pixel so that we can have split up data evenly
    //4 pixels per thread
    //first 16 bytes will contain the number of extra bytes that have been appended to the array
    
    __shared__ unsigned char keys_s[176];
    __shared__ unsigned char sbox_s[256];
    
    sbox_s[threadIdx.x] = sbox_s[threadIdx.x];
    
    if(threadIdx.x < 176){
        keys_s[threadIdx.x] = keys_s[threadIdx.x];
    }
    unsigned char res[8];
    unsigned char mcMatrix[] ={
        0x02, 0x03, 0x01, 0x01,
        0x01, 0x02, 0x03, 0x01,
        0x01, 0x01, 0x02, 0x03,
        0x03, 0x01, 0x01, 0x02
    };
    
    //int tx = (blockIdx.x*BLOCK_SIZE+threadIdx.x)*16;
    int tx = (blockIdx.x*BLOCK_SIZE*threadIdx.x)*32;
    
    unsigned char temp, temp0, temp1;
    
    //make local copy to
    unsigned char l_pic[32];
    
    //printf("!");
    if(tx+32 < size){
    for(int i = 0; i< 32; i++){
        l_pic[i] = pic[tx+i];    
    }
    
    //initial add round key
    for(int i = 0; i < 16; i++){
        temp = keys[i];
        //l_pic[i] = l_pic[i] ^ keys[i];
        l_pic[i] = l_pic[i] ^ temp;
        l_pic[16+i] = l_pic[16+i] ^ temp;
    }
    //The main aes loop
    for(int i = 0; i < 10; i++){
        //sub bytes
        for(int j = 0; j < 32; j++){
            l_pic[j] = sbox_s[l_pic[j]];
        }
        //shift rows
        temp = l_pic[1];
        l_pic[1] = l_pic[1+4];
        l_pic[1+4] = l_pic[1+8];
        l_pic[1+8] = l_pic[1+12];
        l_pic[1+12] = temp;
        
        temp = l_pic[2];
        l_pic[2] = l_pic[2+8];
        l_pic[2+8] = temp;
        temp = l_pic[2+4];
        l_pic[2+4] = l_pic[2+12];
        l_pic[2+12] = temp;
        
        temp = l_pic[3];
        l_pic[3] = l_pic[3+12];
        l_pic[3+12] = l_pic[3+8];
        l_pic[3+8] = l_pic[3+4];
        l_pic[3+4] = temp;
        
        //2nd matrix
        temp = l_pic[16+1];
        l_pic[16+1] = l_pic[16+1+4];
        l_pic[16+1+4] = l_pic[16+1+8];
        l_pic[16+1+8] = l_pic[16+1+12];
        l_pic[16+1+12] = temp;
        
        temp = l_pic[16+2];
        l_pic[16+2] = l_pic[16+2+8];
        l_pic[16+2+8] = temp;
        temp = l_pic[16+2+4];
        l_pic[16+2+4] = l_pic[16+2+12];
        l_pic[16+2+12] = temp;
        
        temp = l_pic[16+3];
        l_pic[16+3] = l_pic[16+3+12];
        l_pic[16+3+12] = l_pic[16+3+8];
        l_pic[16+3+8] = l_pic[16+3+4];
        l_pic[16+3+4] = temp;
        
        //mix columns
        if(i<9){
            for(int j=0; j < 4; j++){
                res[0] = 0x00;
                res[1] = 0x00;
                res[2] = 0x00;
                res[3] = 0x00;
                res[4] = 0x00;
                res[5] = 0x00;
                res[6] = 0x00;
                res[7] = 0x00;
                
                for(int k = 0; k < 4; k++){
                    //res[0] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[0+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1>>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[0] = res[0] ^ temp;
                    
                    temp = 0x00; 
                    temp0 = mcMatrix[0+k];
                    temp1 = l_pic[16+(j*4+k)];
                    for(;temp1;temp1>>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[4+0] = res[4+0] ^ temp;
                    
                    //res[1] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[4+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[1] = res[1] ^ temp;
                    
                    temp = 0x00; 
                    temp0 = mcMatrix[4+k];
                    temp1 = l_pic[16+(j*4+k)];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[4+1] = res[4+1] ^ temp;
                    
                    //res[2] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[8+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[2] = res[2] ^ temp;
                    
                    temp = 0x00; 
                    temp0 = mcMatrix[8+k];
                    temp1 = l_pic[16+(j*4+k)];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[4+2] = res[4+2] ^ temp;
                    
                    //res[3] gf2_mul
                    temp = 0x00; 
                    temp0 = mcMatrix[12+k];
                    temp1 = l_pic[j*4+k];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[3] = res[3] ^ temp;
                    
                    temp = 0x00; 
                    temp0 = mcMatrix[12+k];
                    temp1 = l_pic[16+(j*4+k)];
                    for(;temp1;temp1 >>=1){
                        if(temp1 & 1)
                            temp ^= temp0;
                        if(temp0 & 0x80)
                            temp0 = (temp0 << 1) ^ 0x1b;
                        else
                            temp0 <<= 1;
                    }
                    res[4+3] = res[4+3] ^ temp;
                }
                l_pic[(j*4)+0] = res[0];
                l_pic[(j*4)+1] = res[1];
                l_pic[(j*4)+2] = res[2];
                l_pic[(j*4)+3] = res[3];
                
                l_pic[16+((j*4)+0)] = res[4];
                l_pic[16+((j*4)+1)] = res[5];
                l_pic[16+((j*4)+2)] = res[6];
                l_pic[16+((j*4)+3)] = res[7];
            }
            
        }
        //add round key
        for(int j = 0; j < 16; j++){
            temp = keys_s[((i+1)*16)+j];
            l_pic[j] = l_pic[j] ^ temp;
            l_pic[16+j] = l_pic[16+j] ^ temp;
        }
       //end of of current round of encryption
    }
    //end of all encryption
    //move local bytes back to global memory
    
    for(int i = 0; i < 32; i++){
        r_pic[tx+i] =  l_pic[i];
    }
    }
}



void generateRoundKeys(unsigned char*el_keys, unsigned char *el_sbox, unsigned char *rconTable){
    int currKeyIdx;
    int prevKeyIdx;
    for(int i = 0; i < 10; i++){
        
        currKeyIdx = 16*(i+1);
        prevKeyIdx = 16*i;
        
        el_keys[currKeyIdx+0] = el_keys[prevKeyIdx+13];
        el_keys[currKeyIdx+1] = el_keys[prevKeyIdx+14];
        el_keys[currKeyIdx+2] = el_keys[prevKeyIdx+15];
        el_keys[currKeyIdx+3] = el_keys[prevKeyIdx+12];
        
        el_keys[currKeyIdx+0] = el_sbox[el_keys[currKeyIdx+0]];
        el_keys[currKeyIdx+1] = el_sbox[el_keys[currKeyIdx+1]];
        el_keys[currKeyIdx+2] = el_sbox[el_keys[currKeyIdx+2]];
        el_keys[currKeyIdx+3] = el_sbox[el_keys[currKeyIdx+3]];
        
        el_keys[(16*(i+1))+0] = el_keys[(16*i)+0] ^ el_keys[(16*(i+1))+0] ^ rconTable[i];
        el_keys[(16*(i+1))+1] = el_keys[(16*i)+1] ^ el_keys[(16*(i+1))+1];
        el_keys[(16*(i+1))+2] = el_keys[(16*i)+2] ^ el_keys[(16*(i+1))+2];
        el_keys[(16*(i+1))+3] = el_keys[(16*i)+3] ^ el_keys[(16*(i+1))+3];
        
        for(int j = 1; j < 4; j++){
            el_keys[currKeyIdx+(0+(4*j))] = el_keys[currKeyIdx+(0+(4*(j-1)))] ^ el_keys[prevKeyIdx+(0+(4*j))];
            el_keys[currKeyIdx+(1+(4*j))] = el_keys[currKeyIdx+(1+(4*(j-1)))] ^ el_keys[prevKeyIdx+(1+(4*j))];
            el_keys[currKeyIdx+(2+(4*j))] = el_keys[currKeyIdx+(2+(4*(j-1)))] ^ el_keys[prevKeyIdx+(2+(4*j))];
            el_keys[currKeyIdx+(3+(4*j))] = el_keys[currKeyIdx+(3+(4*(j-1)))] ^ el_keys[prevKeyIdx+(3+(4*j))];
        }
        
    }    
}


void sequentialAES(unsigned char* input, unsigned char* output, unsigned char* keys, unsigned char* sbox, int size){
    
    unsigned char mcMatrix[] = {
        0x02, 0x03, 0x01, 0x01,
        0x01, 0x02, 0x03, 0x01,
        0x01, 0x01, 0x02, 0x03,
        0x03, 0x01, 0x01, 0x02
    };
    
    unsigned char currInput[16];
    unsigned char temp, temp0, temp1;
    unsigned char res[4];
    
    
    int currBaseIdx = 0;
    
    for(int i = 0; i < (size/16) ;i++){
    //for(int i = 0; i < 1; i++){
        currBaseIdx = i*16;
        
        //load in current input section
        for(int j = 0; j < 16; j++){
            currInput[j] = input[currBaseIdx+j];
        }
        /*printf("\ninitial input in seq function: ");
        for(int j=0; j < 16; j++){
            printf("%x ", currInput[j]);
        }*/
        
        //initial add round key
        for(int j = 0; j < 16; j++){
            currInput[j] = currInput[j] ^ keys[j];
        }
        
        /*printf("\nafter init add round key (seq): ");
        for(int j=0; j < 16; j++){
            printf("%x ", currInput[j]);
        }*/
        //begin main aes loop (10 rounds)
        for(int j = 0; j < 10; j++){
            
            //sub bytes
            for(int k = 0; k < 16; k++){
                currInput[k] = sbox[currInput[k]];
            }
            
            /*printf("\nafter sub bytes (seq): ");
            for(int k=0; k < 16; k++){
                printf("%x ", currInput[k]);
            } */
            //shift rows
            temp = currInput[1];
            currInput[1] = currInput[1+4];
            currInput[1+4] = currInput[1+8];
            currInput[1+8] = currInput[1+12];
            currInput[1+12] = temp;
            
            temp = currInput[2];
            currInput[2] = currInput[2+8];
            currInput[2+8] = temp;
            temp = currInput[2+4];
            currInput[2+4] = currInput[2+12];
            currInput[2+12] = temp;
            
            temp = currInput[3];
            currInput[3] = currInput[3+12];
            currInput[3+12] = currInput[3+8];
            currInput[3+8] = currInput[3+4];
            currInput[3+4] = temp;
            
            /*printf("\nafter shift rows (seq): ");
            for(int k=0; k < 16; k++){
                printf("%x ", currInput[k]);
            }*/
            //mix columns
            
            if(j<9){
                
                for(int k = 0; k < 4; k++){
                    
                    res[0] = 0x00;
                    res[1] = 0x00;
                    res[2] = 0x00;
                    res[3] = 0x00;
                    
                    for(int l = 0; l < 4; l++){
                    
                        temp = 0x00; 
                        temp0 = mcMatrix[0+l];
                        temp1 = currInput[k*4+l];
                        
                        for(;temp1;temp1>>=1){
                            if(temp1 & 1)
                                temp ^= temp0;
                            if(temp0 & 0x80)
                                temp0 = (temp0 << 1) ^ 0x1b;
                            else
                                temp0 <<= 1;
                        }
                        res[0] = res[0] ^ temp;
                    
                    //res[1] gf2_mul
                        temp = 0x00; 
                        temp0 = mcMatrix[4+l];
                        temp1 = currInput[k*4+l];
                        
                        for(;temp1;temp1 >>=1){
                            if(temp1 & 1)
                                temp ^= temp0;
                            if(temp0 & 0x80)
                                temp0 = (temp0 << 1) ^ 0x1b;
                            else
                                temp0 <<= 1;
                        }
                        res[1] = res[1] ^ temp;
                    
                    //res[2] gf2_mul
                        temp = 0x00; 
                        temp0 = mcMatrix[8+l];
                        temp1 = currInput[k*4+l];
                        
                        for(;temp1;temp1 >>=1){
                            if(temp1 & 1)
                                temp ^= temp0;
                            if(temp0 & 0x80)
                                temp0 = (temp0 << 1) ^ 0x1b;
                            else
                                temp0 <<= 1;
                        }
                        res[2] = res[2] ^ temp;
                    
                    //res[3] gf2_mul
                        temp = 0x00; 
                        temp0 = mcMatrix[12+l];
                        temp1 = currInput[k*4+l];
                        
                        for(;temp1;temp1 >>=1){
                            if(temp1 & 1)
                                temp ^= temp0;
                            if(temp0 & 0x80)
                                temp0 = (temp0 << 1) ^ 0x1b;
                            else
                                temp0 <<= 1;
                        }
                        res[3] = res[3] ^ temp;
                    }
                    
                    currInput[(k*4)+0] = res[0];
                    currInput[(k*4)+1] = res[1];
                    currInput[(k*4)+2] = res[2];
                    currInput[(k*4)+3] = res[3];
                }
            }
            
            /*printf("\nafter mis cols (seq): ");
            for(int k=0; k < 16; k++){
                printf("%x ", currInput[k]);
            }*/
            
            //add round key
            for(int k = 0; k < 16; k++){
                currInput[k] = currInput[k] ^ keys[(j+1)*16+k];
            }
            
            /*printf("\nafter add round key (seq): ");
            for(int k=0; k < 16; k++){
                printf("%x ", currInput[k]);
            }*/
            
        }
        //move final encrypted 16 byte block t correcto output index
        for(int j = 0; j < 16; j++){
            output[currBaseIdx+j] = currInput[j];
        }
    }
}

void generateRoundKey(unsigned char* prevKey, unsigned char* newKey, int round, unsigned char* sbox,
    unsigned char* rconTable){
    
    //rotWord    
    newKey[0] = prevKey[13];
    newKey[1] = prevKey[14];
    newKey[2] = prevKey[15];
    newKey[3] = prevKey[12];
    
    //subword
    newKey[0] = sbox[newKey[0]];
    newKey[1] = sbox[newKey[1]];
    newKey[2] = sbox[newKey[2]];
    newKey[3] = sbox[newKey[3]];
    
    //rcon
    newKey[0] = prevKey[0] ^ newKey[0] ^ rconTable[round] ;
    newKey[1] = prevKey[1] ^ newKey[1];
    newKey[2] = prevKey[2] ^ newKey[2];
    newKey[3] = prevKey[3] ^ newKey[3];
    
    for(int i = 1; i< 4; i++){
        newKey[0+(4*i)] = newKey[0+(4*(i-1))] ^ prevKey[0+(4*i)];
        newKey[1+(4*i)] = newKey[1+(4*(i-1))] ^ prevKey[1+(4*i)];
        newKey[2+(4*i)] = newKey[2+(4*(i-1))] ^ prevKey[2+(4*i)];
        newKey[3+(4*i)] = newKey[3+(4*(i-1))] ^ prevKey[3+(4*i)];
    }
    
}

int main (int argc, char **argv){
    unsigned char sbox[] ={
        0x63, 0x7c, 0x77, 0x7b, 0xf2, 0x6b, 0x6f, 0xc5, 0x30, 0x01, 0x67, 0x2b, 0xfe, 0xd7, 0xab, 0x76,
        0xca, 0x82, 0xc9, 0x7d, 0xfa, 0x59, 0x47, 0xf0, 0xad, 0xd4, 0xa2, 0xaf, 0x9c, 0xa4, 0x72, 0xc0,
        0xb7, 0xfd, 0x93, 0x26, 0x36, 0x3f, 0xf7, 0xcc, 0x34, 0xa5, 0xe5, 0xf1, 0x71, 0xd8, 0x31, 0x15,
        0x04, 0xc7, 0x23, 0xc3, 0x18, 0x96, 0x05, 0x9a, 0x07, 0x12, 0x80, 0xe2, 0xeb, 0x27, 0xb2, 0x75,
        0x09, 0x83, 0x2c, 0x1a, 0x1b, 0x6e, 0x5a, 0xa0, 0x52, 0x3b, 0xd6, 0xb3, 0x29, 0xe3, 0x2f, 0x84,
        0x53, 0xd1, 0x00, 0xed, 0x20, 0xfc, 0xb1, 0x5b, 0x6a, 0xcb, 0xbe, 0x39, 0x4a, 0x4c, 0x58, 0xcf,
        0xd0, 0xef, 0xaa, 0xfb, 0x43, 0x4d, 0x33, 0x85, 0x45, 0xf9, 0x02, 0x7f, 0x50, 0x3c, 0x9f, 0xa8,
        0x51, 0xa3, 0x40, 0x8f, 0x92, 0x9d, 0x38, 0xf5, 0xbc, 0xb6, 0xda, 0x21, 0x10, 0xff, 0xf3, 0xd2,
        0xcd, 0x0c, 0x13, 0xec, 0x5f, 0x97, 0x44, 0x17, 0xc4, 0xa7, 0x7e, 0x3d, 0x64, 0x5d, 0x19, 0x73,
        0x60, 0x81, 0x4f, 0xdc, 0x22, 0x2a, 0x90, 0x88, 0x46, 0xee, 0xb8, 0x14, 0xde, 0x5e, 0x0b, 0xdb,
        0xe0, 0x32, 0x3a, 0x0a, 0x49, 0x06, 0x24, 0x5c, 0xc2, 0xd3, 0xac, 0x62, 0x91, 0x95, 0xe4, 0x79,
        0xe7, 0xc8, 0x37, 0x6d, 0x8d, 0xd5, 0x4e, 0xa9, 0x6c, 0x56, 0xf4, 0xea, 0x65, 0x7a, 0xae, 0x08,
        0xba, 0x78, 0x25, 0x2e, 0x1c, 0xa6, 0xb4, 0xc6, 0xe8, 0xdd, 0x74, 0x1f, 0x4b, 0xbd, 0x8b, 0x8a,
        0x70, 0x3e, 0xb5, 0x66, 0x48, 0x03, 0xf6, 0x0e, 0x61, 0x35, 0x57, 0xb9, 0x86, 0xc1, 0x1d, 0x9e,
        0xe1, 0xf8, 0x98, 0x11, 0x69, 0xd9, 0x8e, 0x94, 0x9b, 0x1e, 0x87, 0xe9, 0xce, 0x55, 0x28, 0xdf,
        0x8c, 0xa1, 0x89, 0x0d, 0xbf, 0xe6, 0x42, 0x68, 0x41, 0x99, 0x2d, 0x0f, 0xb0, 0x54, 0xbb, 0x16
    };
    
    unsigned char rconTable[]={0x01, 0x02, 0x04, 0x08, 0x10, 0x20, 0x40, 0x80, 0x1b, 0x36};
    
    
    int thresh = 8000;
    char *filename;
    filename = strdup(DEFAULT_FILENAME);
    
    if(argc > 1){
        if(argc == 2){
            filename = strdup(argv[1]);
        } 
        fprintf(stderr, "file %s    threshold %d\n", filename, thresh);
    }
    
    int xsize, ysize, maxval;
    unsigned char *pic = read_ppm(filename, &xsize, &ysize, &maxval);
    
    printf("xsize: %d, ysize: %d, maxval: %d\n", xsize, ysize, maxval);
    //zero out the pic first because you need to know how much of the array is actually being used
    
    int numbytes = xsize * ysize * 3 * sizeof(int);
    //int numbytesp = xsize *ysize * 3 *sizeof(int)*3;
    
    //int *result = (int *) malloc(numbytes);
    
    //unsigned char *resultc = (unsigned char *) malloc(numbytes*3);
    
    printf("size of int: %ld", sizeof(int));
    unsigned char *resultc = (unsigned char *) malloc(xsize*ysize*3);
    if(!resultc){
        fprintf(stderr, "sobel() unable to malloc %d bytes\n", numbytes); 
        exit(-1);
    }
    
    //memset(resultc, 0x00, numbytes*3);
    memset(resultc, 0x00, xsize*ysize*3);
    
    //memset(resultc, 0, xsize * ysize*sizeof(int)*3);
    
    printf("got here yurd baby 3\n");
    //int i, j, offset;
    //char temp[3];
    //please kill me
    printf("got here 1 14/16\n");
    printf("sizes: %d, %d\n",xsize,ysize);
    
    
    /*for(i = 0; i< ysize; i++){
        for(j = 0; j < xsize; j++){
            offset = i * xsize + j;
            //offset = offset +1 -1;
            
            //printf("%x,%x,%x ", pic[offset*3],pic[offset*3+1],pic[offset*3+2]);
            resultc[offset*3] = pic[offset*3];
            //resultc[offset*3] = 0x00;
            resultc[offset*3+1] = pic[offset*3+1];
            //resultc[offset*3+1] = 0x00;
            //resultc[offset*3+2] = pic[offset*3+2];
            resultc[offset*3+2] = pic[offset*3+2];
            resultc[offset*3+2] = 0x80;
            //result = (int*)(pic[offset*3])
            //((*int)&temp[0]) = (pic[offset]);
            //printf("%x,%x ",temp[0],temp[1]);
            //itoa(pic[offset], temp, 16);
            //printf("%s, ",temp);
            
            //temp[0]= *((char*)pic[offset]);
            //temp[1]= *(((char*)pic[offset+1]));
            //temp[2]= *((char*)pic[offset]);
            //printf("(%x,%x,%x) ",temp[0], temp[1], temp[2]);
            
            
                
            //magnitude = sum1*sum1 + sum2*sum2;
            
            //if(magnitude > thresh)
                //result[offset] = 128;
            //else
                //result[offset] = 0;
        }
        //printf("\n");
    }*/
    /*printf("\nxsize * ysize * 3 * sizeof(int)*3 - 1: %x\n xsize * ysize * 3 * sizeof(int)-1: %x\n xsize * ysize*3 -1:%x\n(xsize * ysize) -1: %x\n",
        resultc[(xsize*ysize*3*sizeof(int)*3)-5], resultc[(xsize*ysize*3*sizeof(int))-5], 
        resultc[(xsize*ysize*3)-5], resultc[(xsize*ysize)-5]);*/
        
       //////////   BEGINNING OF SEQUENTIAL AES    ////////////////////////////////////////////////////////////////// 
        
        
        printf("beginning (fixed?): ");
        for(int i = 0; i < 16; i++){
            printf("%x ", pic[i]);
        }
        printf("\n");
    //GENERATE KEY!
    
    //unsigned char roundKeys[11][16];
    unsigned char roundKeys2 [176]; //11*16 = 176 quick maths
    //unsigned char roundKeyHolder...
    
    /*unsigned char testKey[] = {
        0x21, 0x8f, 0x0d, 0x3e, 
        0x95, 0x3c, 0xb0, 0x19,
        0xd1, 0x63, 0x4b, 0xf4,
        0x12, 0xab, 0x8f, 0x3c
    };*/
    
    unsigned char testKey2[] = {
        0x2b, 0x7e, 0x15, 0x16,
        0x28, 0xae, 0xd2, 0xa6,
        0xab, 0xf7, 0x15, 0x88,
        0x09, 0xcf, 0x4f, 0x3c
    };
   /* 
    for(int i = 0; i< 16; i++){
        roundKeys[0][i] = testKey2[i];
    }*/
    for(int i = 0; i < 16; i++){
        roundKeys2[i] = testKey2[i];
    }
    
    generateRoundKeys(roundKeys2, sbox, rconTable);
    
    for(int i = 0; i < 11; i++){
        printf("k%d: ",i);
        for(int j = 0; j < 16; j++){
            printf("%x ",roundKeys2[(i*16)+j]);
        }
        printf("\n");
    }
    
    //clock_t start_s, end_s;
    //double cpu_time_used;
    
    //start_s=clock();
    cudaEvent_t start_s, stop_s;
    float deltaTime_s;
    
    cudaEventCreate(&start_s);
    cudaEventCreate(&stop_s);
    cudaEventRecord(start_s, 0);
    
    sequentialAES(pic, resultc, roundKeys2, sbox, xsize*ysize*3);
    
    cudaEventRecord(stop_s,0);
    cudaEventSynchronize(stop_s);
    cudaEventElapsedTime(&deltaTime_s,start_s,stop_s);
    cudaEventDestroy(start_s);
    cudaEventDestroy(start_s);
    
    //end_s= clock();
    //cpu_time_used = ((double) (end_s-start_s)) /CLOCKS_PER_SEC;
    
    //printf("Sequential AES time = %f\n", cpu_time_used);
    printf("Delta Time (seq) = %f\n", deltaTime_s);
    
    printf("sequential first 16 bytes result: ");
    for(int i =0; i < 16; i++){
        printf("%x ", resultc[i]);
    }
    printf("\n");
    //first round of generate new key
    //generateRoundKey(testKey2, roundKeys[0], 0, sbox, rconTable);
    /*for(int i = 1; i < 11; i++){
        generateRoundKey(roundKeys[i-1], roundKeys[i], i-1,sbox,rconTable);
    }*/
    //printing roundkeys
    /*for(int i = 0; i < 11; i++){
        for(int j = 0; j < 16; j++){
            printf("%x ", roundKeys[i][j]);
        }
        printf("\n");
    }*/
    
    
    unsigned char* h_output = (unsigned char*)malloc(xsize*ysize*3);
    
    unsigned char* test_host_input = (unsigned char*)malloc(16*sizeof(unsigned char));
    unsigned char* test_host_output = (unsigned char*)malloc(16*sizeof(unsigned char));
    unsigned char* h_sbox = (unsigned char*)malloc(16*16*sizeof(unsigned char));
    unsigned char* h_roundKeys = (unsigned char*)malloc(176*sizeof(unsigned char));
    //unsigned char* h_rconTable = (unsigned char*)malloc(10*sizeof(unsigned char));
    
    unsigned char* d_input;
    unsigned char* d_output;
    
    unsigned char* test_device_input;
    unsigned char* test_device_output;
    unsigned char* d_sbox;
    unsigned char* d_roundKeys;
    //unsigned char* d_rconTable;
    
    cudaMalloc(&d_input, xsize*ysize*3);
    cudaMalloc(&d_output, xsize*ysize*3);
    
    cudaMalloc(&test_device_input, 16*sizeof(unsigned char));
    cudaMalloc(&test_device_output, 16*sizeof(unsigned char));
    cudaMalloc(&d_sbox, 16*16*sizeof(unsigned char));
    cudaMalloc(&d_roundKeys, 176*sizeof(unsigned char));
    //cudaMalloc(&d_rconTable, 10*sizeof(unsigned char));
    
    //initialize test input
    
    test_host_input[0] = 0x56;
    test_host_input[1] = 0xf6;
    test_host_input[2] = 0x4b;
    test_host_input[3] = 0xb5;
    test_host_input[4] = 0x9c;
    test_host_input[5] = 0x7d;
    test_host_input[6] = 0x17;
    test_host_input[7] = 0xa4;
    test_host_input[8] = 0x9b;
    test_host_input[9] = 0x08;
    test_host_input[10] = 0x90;
    test_host_input[11] = 0x02;
    test_host_input[12] = 0x50;
    test_host_input[13] = 0xeb;
    test_host_input[14] = 0xbb;
    test_host_input[15] = 0x3f;
    
    for(int i = 0; i < 16*16; i++){
        h_sbox[i] = sbox[i];
    }
    for(int i = 0; i < 176; i++){
        h_roundKeys[i] = roundKeys2[i];
    }
    /*for(int i = 0; i < 10; i++){
        h_rconTable[i] = rconTable[i];
    }*/
    
    cudaMemcpy(d_input, pic, xsize*ysize*3, cudaMemcpyHostToDevice);
    
    
    cudaMemcpy(test_device_input, test_host_input,16*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_sbox, h_sbox,16*16*sizeof(unsigned char), cudaMemcpyHostToDevice);
    cudaMemcpy(d_roundKeys, h_roundKeys, 176*sizeof(unsigned char),cudaMemcpyHostToDevice);
    //cudaMemcpy(d_rconTable, h_rconTable, cudaMemcpyHostToDevice);
    
    cudaEvent_t start, stop;
    float deltaTime;
    
    cudaEventCreate(&start);
    cudaEventCreate(&stop);
    cudaEventRecord(start, 0);
    
    aes_moreWork<<<(int)floor((float)(xsize*ysize*3)/(BLOCK_SIZE*16)),BLOCK_SIZE>>>(d_input, d_output, d_roundKeys, d_sbox, xsize*ysize*3);
    
    //aes_moreWork<<<(int)floor((float)(xsize*ysize*3)/(BLOCK_SIZE*32)),BLOCK_SIZE>>>(d_input, d_output, d_roundKeys, d_sbox, xsize*ysize*3);
    printf("\n");
    //aes_moreWork<<<1,1>>>(d_input, d_output, d_roundKeys, d_sbox);
    
    cudaEventRecord(stop,0);
    cudaEventSynchronize(stop);
    cudaEventElapsedTime(&deltaTime,start,stop);
    cudaEventDestroy(start);
    cudaEventDestroy(start);
    //aes<<<1,1>>>(test_device_input,test_device_output, d_roundKeys, d_sbox);
    
    cudaMemcpy(h_output,d_output, xsize*ysize*3, cudaMemcpyDeviceToHost);
    //cudaMemcpy(test_host_output, test_device_output, 16*sizeof(unsigned char), cudaMemcpyDeviceToHost);
    
    printf("Delta Time (more work) = %f\n", deltaTime);
    
    printf("\nFinal Output!: ");
    for(int i = 0; i < 16; i++){
        printf("%x ",h_output[i]);
    }
    
    printf("got here 2 son\n");
    write_ppm((unsigned char *) "result8sendhelp.ppm", xsize, ysize, 255, resultc);
    write_ppm((unsigned char *) "encrypted_gary.ppm", xsize,ysize, 255, h_output);
    
    int mismatchCounter = 0;
    for(int i =0 ; i < xsize*ysize*3; i++){
        if(resultc[i] != h_output[i]){
            mismatchCounter++;
        }
    }
    printf("mismatchCounter: %d\n", mismatchCounter);
    
    free(pic);
    free(resultc);
    free(test_host_input);
    free(test_host_output);
    free(h_sbox);
    free(h_output);
    //free(h_)
    
    cudaFree(test_device_input);
    cudaFree(test_device_output);
    cudaFree(d_sbox);
    cudaFree(d_output);
    cudaFree(d_input);
    
    fprintf(stderr, "aes done\n");
}


