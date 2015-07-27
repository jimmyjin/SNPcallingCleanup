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
import com.snp.callings.tools.VCFLookUpDataSource;

// The main class to handle the cleanup
public class JackyResultCleanup {

	private String snpVcf;
	private String snpCsv;
	private String indelsVcf;
	private String indelsCsv;
	private String finalOutput;
	private String combinedTempFile;

	
	public JackyResultCleanup(String subjectId, String directory)
	{
		snpVcf = directory + "/" + subjectId + "_snps_annovar.hg19_multianno.vcf";
		snpCsv = directory + "/" + subjectId + "_snps_annovar.hg19_multianno.csv";
		indelsVcf = directory + "/" + subjectId + "_indels_annovar.hg19_multianno.vcf";
		indelsCsv = directory + "/" + subjectId + "_indels_annovar.hg19_multianno.csv";
		finalOutput = directory + "/" + subjectId + "_variant_calling_final_result.csv";
		combinedTempFile = directory + "/" + subjectId + "_combined_temp.csv";
	}

	public void appendColumns()
	{
		try {
			CsvReader snpCsvReader = new CsvReader(snpCsv);
			VCFLookUpDataSource snpVcfDataStore = new VCFLookUpDataSource(snpVcf);
			CsvReader indelsCsvReader = new CsvReader(indelsCsv);
			VCFLookUpDataSource indelsVcfDataStore = new VCFLookUpDataSource(indelsVcf);
			CsvWriter combinedOutput =new CsvWriter(combinedTempFile,',',Charset.forName("UTF-8"));
			
			// read csv header
			snpCsvReader.readHeaders();
			indelsCsvReader.readHeaders();
			List<String> headers = new ArrayList<String>();
			Collections.addAll(headers, snpCsvReader.getHeaders());
			headers.add("zygosity");
			headers.add("coverage");
			headers.add("src"); // for debug purpose, to rename todo
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
			CsvReader combinedCsvReader = new CsvReader(combinedTempFile);
			CsvWriter finalOutputWriter =new CsvWriter(finalOutput,',',Charset.forName("UTF-8"));
			combinedCsvReader.readHeaders();
			finalOutputWriter.writeRecord(header);
			
			while (combinedCsvReader.readRecord())
			{
				String chr = combinedCsvReader.get("Chr");
				if (chr.startsWith("chr")){
					chr = chr.substring(3);
				}
				String[] content = {chr + ":" + clean(combinedCsvReader.get("Start")), clean(combinedCsvReader.get("snp138")), clean(combinedCsvReader.get("Gene.refGene")), 
						clean(combinedCsvReader.get("coverage")), "",  clean(combinedCsvReader.get("GeneDetail.refGene")), clean(combinedCsvReader.get("AAChange.refGene")), clean(combinedCsvReader.get("ExonicFunc.refGene")), "",
						clean(combinedCsvReader.get("zygosity")), clean(combinedCsvReader.get("SIFT_pred")), clean(combinedCsvReader.get("Polyphen2_HDIV_pred")), clean(combinedCsvReader.get("Polyphen2_HVAR_pred")),clean(combinedCsvReader.get("Ref")),
						clean(combinedCsvReader.get("Alt")), clean(combinedCsvReader.get("Func.refGene"))};
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

	private void generateOutput(CsvReader CsvReader, CsvWriter combinedOutput, VCFLookUpDataSource vcfDataStore) throws IOException
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
				if (af.startsWith("1.00"))
				{
					zygosity = "homozygous";
				} else if (af.startsWith("0.500")) {
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
	
	private String clean(String input)
	{
		if (input.equals("."))
		{
			return "";
		}
		return input.replace(',', ';');
	}

	public static void main( String[] args )
	{
//		if (args.length != 3)
//		{
//			System.out.println("The use of the tool");
//			System.exit(1);
//		}
		JackyResultCleanup j = new JackyResultCleanup("WGC033595U_combined", "c:/Users/sjin1/Downloads/jackyfiles/jackyfiles");
		j.appendColumns();
		j.mappingColumns();
	}

}
