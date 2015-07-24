/**
 * The Clean-up tool to get final result for jacky's SNPcalling tool
 *
 * @Author Shengmin (Jimmy) Jin
 * Created on 7/24/2015.
 */

package com.snp.callings.CleanUp;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.ArrayList;
import java.util.Collections;
import java.util.List;
import java.util.Map;

import com.csvreader.CsvReader;
import com.csvreader.CsvWriter;
import com.snp.callings.tools.VCFDataStore;


public class JackyResultCleanup {

	private String subjectId;
	private String directory;
	private String snpVcf;
	private String snpCsv;
	private String indelsVcf;
	private String indelsCsv;
	private String finalOutput;
	private String combinedFile;

	
	public JackyResultCleanup(String subjectId, String directory)
	{
		this.subjectId = subjectId;
		this.directory = directory;
		snpVcf = directory + "/" + subjectId + "_snps_annovar.hg19_multianno.vcf";
		snpCsv = directory + "/" + subjectId + "_snps_annovar.hg19_multianno.csv";
		indelsVcf = directory + "/" + subjectId + "_indels_annovar.hg19_multianno.vcf";
		indelsCsv = directory + "/" + subjectId + "_indels_annovar.hg19_multianno.csv";
		finalOutput = directory + "/" + subjectId + "_variant_calling_final_result.csv";
		combinedFile = directory + "/" + subjectId + "_combined_temp.csv";
	}

	public void appendColumns()
	{

		// snps
		try {
			CsvReader snpCsvReader = new CsvReader(snpCsv);
			VCFDataStore snpVcfDataStore = new VCFDataStore(snpVcf);
			CsvReader indelsCsvReader = new CsvReader(indelsCsv);
			VCFDataStore indelsVcfDataStore = new VCFDataStore(indelsVcf);
			CsvWriter combinedOutput =new CsvWriter(combinedFile,',',Charset.forName("UTF-8"));
			
			// read csv header
			snpCsvReader.readHeaders();
			indelsCsvReader.readHeaders();
			List<String> headers = new ArrayList<String>();
			Collections.addAll(headers, snpCsvReader.getHeaders());
			headers.add("zygosity");
			headers.add("coverage");
			headers.add("src"); // for debug
			String[] headerstring = new String[headers.size()];
			combinedOutput.writeRecord(headers.toArray(headerstring));

			generateOutput(snpCsvReader, combinedOutput, snpVcfDataStore);
			generateOutput(indelsCsvReader, combinedOutput, indelsVcfDataStore);
			
			snpCsvReader.close();
			indelsCsvReader.close();
			combinedOutput.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}
	
	public void mappingColumns()
	{
		String[] header = {"chr_position","dbsnp_id","gene","coverage","exon_number","hgvs_coding","hgvs_protein","variant_result_type",
				"transcript_accession","zygosity","sift","polyphen2_hdiv","polyphen2_hvar","ref","alt","variant_location"};
		try {
			CsvReader combinedCsvReader = new CsvReader(combinedFile);
			CsvWriter finalOutputWriter =new CsvWriter(finalOutput,',',Charset.forName("UTF-8"));
			combinedCsvReader.readHeaders();
			finalOutputWriter.writeRecord(header);
			
			while (combinedCsvReader.readRecord())
			{
				String chr = combinedCsvReader.get("Chr");
				if (chr.startsWith("chr")){
					chr = chr.substring(3);
				}
				String[] content = {chr + ":" + combinedCsvReader.get("Start"), combinedCsvReader.get("snp138"), combinedCsvReader.get("Gene.refGene"), 
						combinedCsvReader.get("coverage"), "",  combinedCsvReader.get("GeneDetail.refGene"), combinedCsvReader.get("AAChange.refGene"), combinedCsvReader.get("ExonicFunc.refGene"), "",
						combinedCsvReader.get("zygosity"), combinedCsvReader.get("SIFT_pred"), combinedCsvReader.get("Polyphen2_HDIV_pred"), combinedCsvReader.get("Polyphen2_HVAR_pred"),combinedCsvReader.get("Ref"),
						combinedCsvReader.get("Alt"), combinedCsvReader.get("Func.refGene")};
				finalOutputWriter.writeRecord(content);
			}
			
			combinedCsvReader.close();
			finalOutputWriter.close();
		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private void generateOutput(CsvReader CsvReader, CsvWriter combinedOutput, VCFDataStore vcfDataStore) throws IOException
	{
		Map<String, String> info;
		while (CsvReader.readRecord())
		{
			String ID, mapID, zygosity = "", coverage = "", src = "";
			ID = CsvReader.get("snp138");
			mapID = CsvReader.get("Chr") + ":" + CsvReader.get("Start") + ":" + CsvReader.get("Ref") + ":" + 
					CsvReader.get("Alt");
			if (!ID.equals("."))
			{
				info = vcfDataStore.firstLookup(ID);
				src = "1";
				if (info == null)
				{
					info = vcfDataStore.secondaryLookup(mapID);
					src = "2";
				}
			} else {
				info = vcfDataStore.secondaryLookup(mapID);
				src = "2";
			}

			if (info != null)
			{
				String af = info.get("AF");
				if (af.equals("1.00"))
				{
					zygosity = "homozygous";
				} else if (af.equals("0.500")) {
					zygosity = "heterozygous";
				}
				coverage= info.get("DP");
			}
			List<String> contents = new ArrayList<String>();
			Collections.addAll(contents, CsvReader.getValues());
			contents.add(zygosity);
			contents.add(coverage);
			contents.add(src);
			String[] strings = new String[contents.size()];
			combinedOutput.writeRecord(contents.toArray(strings));
		}
	}


	public static void main( String[] args )
	{
		JackyResultCleanup j = new JackyResultCleanup("WGC033595U_combined", "c:/Users/sjin1/Downloads/jackyfiles/jackyfiles");
		j.appendColumns();
		j.mappingColumns();
	}

}
