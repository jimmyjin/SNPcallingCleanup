package com.snp.callings.tools;

import java.io.FileNotFoundException;
import java.io.IOException;
import java.nio.charset.Charset;
import java.util.HashMap;
import java.util.Map;

import com.csvreader.CsvReader;
import com.csvreader.CsvWriter;

public class VCFLookUpDataSource {

	private Map<String, Map<String, String>> primaryMap = new HashMap<String, Map<String, String>>();
	private Map<String, Map<String, String>> secondaryMap = new HashMap<String, Map<String, String>>();
	private Map<String, Integer> header = new HashMap<String, Integer>();

	public VCFLookUpDataSource(String vcfFile)
	{
		generateFromVcf(vcfFile);
	}

	private void  generateFromVcf(String vcfFile){
		try {

			CsvReader reader = new CsvReader(vcfFile, '\t');

			while (reader.readRecord())
			{
				String rec = reader.getRawRecord();
				if (rec.startsWith("##")) continue;
				if (rec.startsWith("#"))
				{
					String[] head_list = rec.substring(1).split("\t");
					reader.setHeaders(head_list);
					for (int i = 0; i < head_list.length; ++i)
					{
						header.put(head_list[i], i);
					}
					continue;
				}

				//String id = reader.get("ID");
				Map<String, String> info = parseInfo(reader.get("INFO"));
				String primary_id = info.get("snp138");
				if (primary_id != "." && primary_id != "")
				{
					primaryMap.put(primary_id, info);
				}
				
				String secondaryID = reader.get("CHROM") + ":" + reader.get("POS") + ":" + reader.get("REF") + ":" + reader.get("ALT");
				secondaryMap.put(secondaryID, parseInfo(reader.get("INFO")));
			}

			reader.close();

		} catch (FileNotFoundException e) {
			e.printStackTrace();
		} catch (IOException e) {
			e.printStackTrace();
		}
	}

	private Map<String, String> parseInfo(String info)
	{
		Map<String, String> m = new HashMap<String, String>();
		String[] attributes = info.split(";");
		for (String a : attributes)
		{
			String[] kvPair = a.split("=");
			if (kvPair.length == 2)
			{
				m.put(kvPair[0], kvPair[1]);
			}
		}
		return m;
	}

	public Map<String, String> firstLookup(String id)
	{
		if (primaryMap.containsKey(id))
		{
			return primaryMap.get(id);
		} else {
			return null;
		}
	}

	public Map<String, String> secondaryLookup(String mapId)
	{
		if (secondaryMap.containsKey(mapId))
		{
			return secondaryMap.get(mapId);
		} else {
			return null;
		}
	}


	public static void main( String[] args )
	{
		VCFLookUpDataSource j = new VCFLookUpDataSource("c:/Users/sjin1/Downloads/jackyfiles/jackyfiles/WGC033595U_combined_snps_annovar.hg19_multianno.vcf");

		System.out.println("primary map");
		for (Map.Entry<String, String> entry : j.firstLookup("rs75062661").entrySet())
		{
			System.out.println(entry.getKey() + "::" + entry.getValue());
		}

	}
}
