#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include "curl.h"
#include "oxtypes.h"
#include "oxexport.h"
#include "jdmath.h"

#define OXARGS (OxVALUE *rtn, OxVALUE *pv, int cArg)

void OXCALL FnThrees OXARGS
{
    int i, j, c, r;

    OxLibCheckType(OX_INT, pv, 0, 1);

    r = OxInt(pv, 0);
    c = OxInt(pv, 1);
    OxLibValMatMalloc(rtn, r, c);

    for (i = 0; i < r; i++)
        for (j = 0; j < c; j++)
            OxMat(rtn, 0)[i][j] = 3;
}


static size_t write_data(void *ptr, size_t size, size_t nmemb, void *stream) {
  size_t written = fwrite(ptr, size, nmemb, (FILE *)stream);
  return written;
  }


void OXCALL FnGet OXARGS {
    OxLibCheckType(OX_STRING,pv,0,1);
 CURL *curl_handle;
      FILE *pagefile;

    curl_global_init(CURL_GLOBAL_ALL);
    curl_handle = curl_easy_init();
    curl_easy_setopt(curl_handle, CURLOPT_URL, OxStr(pv,0));

  curl_easy_setopt(curl_handle, CURLOPT_NOPROGRESS, 1L);
  curl_easy_setopt(curl_handle, CURLOPT_WRITEFUNCTION, write_data);

  pagefile = fopen(OxStr(pv,1), "wb");
  if(pagefile) {
    curl_easy_setopt(curl_handle, CURLOPT_WRITEDATA, pagefile);
    curl_easy_perform(curl_handle);
    fclose(pagefile);
  }
  curl_easy_cleanup(curl_handle);

  }
