#include <jni.h>
#include <stdlib.h>
#include <string.h>
#include "fermi-lite/fml.h"

JNIEXPORT jobject JNICALL
Java_org_broadinstitute_hellbender_utils_fermi_FermiLiteAssembler_createAssemblyData( JNIEnv* env, jclass cls, jobject readBuf ) {

	// build input data structure (reads)
	int32_t* pReadsBuf = (*env)->GetDirectBufferAddress(env, readBuf);
	int32_t nSeqs = *pReadsBuf++;
	bseq1_t* pReads = malloc(nSeqs*sizeof(bseq1_t));
	bseq1_t* pReadsEnd = pReads + nSeqs;
	bseq1_t* pRead;
	char* readData = (char*)pReadsBuf;
	for ( pRead = pReads; pRead != pReadsEnd; ++pRead ) {
		int skip = strlen(readData);
		pRead->l_seq = skip++;
		pRead->seq = strdup(readData);
		readData += skip;
		pRead->qual = strdup(readData);
		readData += skip;
	}

	// assemble the reads
	int32_t nUnitigs;
	fml_opt_t opts;
	fml_opt_init(&opts);
	fml_utg_t* pUnitigs = fml_assemble(&opts, nSeqs, pReads, &nUnitigs);
	fml_utg_t* pUnitigsEnd = pUnitigs + nUnitigs;

	// marshal the output data (unitigs and connections)
	size_t arrSize = 2*sizeof(int32_t); // for the array length and pool offset
	size_t totSize = 0;
	fml_utg_t* pUnitig;
	for ( pUnitig = pUnitigs; pUnitig != pUnitigsEnd; ++pUnitig ) {
		totSize += 2*pUnitig->len; // for the contig sequence and per-base coverage bytes
		// each contig has 3 ints (length, #supporting reads, nConnections)
		// and each of the connections requires an fml_ovlp_t to describe the edge
		arrSize += 3*sizeof(int32_t) + (pUnitig->n_ovlp[0]+pUnitig->n_ovlp[1])*sizeof(fml_ovlp_t);
	}
	totSize += arrSize;
	int32_t* pAsmBuf = malloc(totSize);
	int32_t* pTig = pAsmBuf;
	*pTig++ = nUnitigs; // number of unitigs
	*pTig++ = arrSize;  // offset to byte pool
	// copy the per-unitig data
	for ( pUnitig = pUnitigs; pUnitig != pUnitigsEnd; ++pUnitig ) {
		*pTig++ = pUnitig->len;
		*pTig++ = pUnitig->nsr;
		int32_t nConnections = pUnitig->n_ovlp[0] + pUnitig->n_ovlp[1];
		*pTig++ = nConnections;
		size_t ovlpLen = nConnections*sizeof(fml_ovlp_t);
		memcpy(pTig, pUnitig->ovlp, ovlpLen);
		pTig += ovlpLen/sizeof(int32_t);
	}
	// copy the per-unitig byte data (i.e., sequence and coverage)
	char* pBytes = (char*)pTig;
	for ( pUnitig = pUnitigs; pUnitig != pUnitigsEnd; ++pUnitig ) {
		size_t len = pUnitig->len;
		memcpy(pBytes, pUnitig->seq, len);
		pBytes += len;
		memcpy(pBytes, pUnitig->cov, len);
		pBytes += len;
	}

	// clean up assembly data
	fml_utg_destroy(nUnitigs, pUnitigs);

	// return the output to Java
	jobject result = (*env)->NewDirectByteBuffer(env, pAsmBuf, totSize);
	if ( !result ) free(pAsmBuf);
	return result;
}

JNIEXPORT void JNICALL
Java_org_broadinstitute_hellbender_utils_fermi_FermiLiteAssembler_destroyAssemblyData( JNIEnv* env, jclass cls, jobject asmBuf ) {
	free((*env)->GetDirectBufferAddress(env, asmBuf));
}
